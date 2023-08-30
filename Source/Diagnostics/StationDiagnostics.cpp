/* Copyright 2023 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "StationDiagnostics.H"
#include "ComputeDiagFunctors/StationFunctor.H"
//temporary
#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include "ComputeDiagFunctors/StationParticleFunctor.H"
#include "Diagnostics/Diagnostics.H"
#include "Diagnostics/FlushFormats/FlushFormat.H"

#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

//temporary
#include <ablastr/coarsen/sample.H>
#include <ablastr/utils/Communication.H>
#include <ablastr/utils/SignalHandling.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_CoordSys.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_FabFactory.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MakeType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_ParticleContainerBase.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_VisMF.H>
#include <AMReX_Vector.H>
#include <AMReX_VectorIO.H>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

using namespace amrex::literals;

StationDiagnostics::StationDiagnostics (int i, std::string name)
    : Diagnostics(i, std::move(name))
{
    ReadParameters();
}

void
StationDiagnostics::ReadParameters ()
{
#ifdef WARPX_DIM_RZ
    amrex::Abort("StationDiagnostics is not implemented for RZ, yet");
#endif

    BaseReadParameters();
    const amrex::ParmParse pp_diag_name(m_diag_name);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_format == "plotfile", "<diag>.format must be plotfile");

    std::string station_normal_dir;
    pp_diag_name.get("station_normal_dir", station_normal_dir);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        station_normal_dir == "z", "<diag>.station_normal_dir is supported only for z at the moment");
    if (station_normal_dir == "z") {
        m_station_normal = StationNormalDir::z;
    }

    pp_diag_name.get("station_location",m_station_loc);

    // If station recordings are used to restart simulations, then raw fields are needed
    pp_diag_name.query("plot_raw_fields", m_plot_raw_fields);
    pp_diag_name.query("plot_raw_fields", m_plot_raw_fields_guards);

    pp_diag_name.query("buffer_size",m_buffer_size);
    // for now, number of buffers or in this case, z-locations is assumed to be 1
    // This is used to allocate the number of output multi-level multifabs, m_mf_output
    m_num_buffers = 1;

    std::string intervals_string_vec = "0";
    pp_diag_name.query("intervals", intervals_string_vec);
    m_intervals = utils::parser::SliceParser(intervals_string_vec);

    m_varnames = {"Ex", "Ey", "Ez", "Bx", "By", "Bz"};
    m_file_prefix = "diags/" + m_diag_name;
    pp_diag_name.query("file_prefix", m_file_prefix);

    auto& warpx = WarpX::GetInstance();
    MultiParticleContainer& mpc = warpx.GetPartContainer();
    int write_species = 1;
    pp_diag_name.query("write_species", write_species);
    if (m_output_species_names.size() == 0 and write_species == 1)
        m_output_species_names = mpc.GetSpeciesNames();

    bool m_record_particles = false;
    if (m_output_species_names.size() > 0 and write_species == 1) {
        m_record_particles = true;
    }
    if (m_record_particles) {
        // we use the same logic. I dont think we need a separate variable. but might be more explicit
        mpc.SetDoBackTransformedParticles(m_record_particles);
        for (auto const& species : m_output_species_names) {
            mpc.SetDoBackTransformedParticles(species, m_record_particles);
        }
    }
    const int num_stationdiag_buffers = 1;
    m_particles_buffer.resize(num_stationdiag_buffers);
    m_totalParticles_in_buffer.resize(num_stationdiag_buffers);
}

bool
StationDiagnostics::DoDump (int step, int /* i_buffer*/, bool force_flush)
{
    // Determine criterion to dump output
    return ( (m_slice_counter == m_buffer_size) || force_flush || m_last_timeslice_filled);
}

bool
StationDiagnostics::DoComputeAndPack (int step, bool force_flush)
{
    // we may have to compute and pack everytimestep, but only store at user-defined intervals
    return ( step>=0 );
}

void
StationDiagnostics::InitializeFieldFunctors (int lev)
{
    // need to define a functor that will take a slice of the fields and store it in a given location
    // in the multifab based on current time
    // In this function, we will define m_all_field_functors based on that functor
    auto & warpx = WarpX::GetInstance();
    const int num_stationdiag_functors = 1;
    m_all_field_functors[lev].resize(num_stationdiag_functors);
    int ncomp = 6;
    for (int i = 0; i < num_stationdiag_functors; ++i)
    {
        m_all_field_functors[lev][i] = std::make_unique<StationFunctor>(
                                           warpx.get_array_EBfield_fp(lev), m_station_loc, lev,
                                           m_crse_ratio, ncomp
                                       );
    }
}

void
StationDiagnostics::InitializeBufferData (int i_buffer, int lev, bool restart)
{
    // Define boxArray, dmap, and initialize output multifab
    // This will have the extension in x-y at a given location z, and the third dimension will be time

    auto & warpx = WarpX::GetInstance();
    // in the station normal direction, with time, we need to determine number of points
    // hi : (t_max - t_min )/dt
    amrex::Box domain = (warpx.boxArray(lev)).minimalBox();
    domain.setSmall(WARPX_ZINDEX, 0);
    domain.setBig(WARPX_ZINDEX, (m_buffer_size - 1));
    m_buffer_box = domain;
    amrex::BoxArray diag_ba;
    diag_ba.define(m_buffer_box);
    amrex::BoxArray ba = diag_ba.maxSize(256);
    amrex::DistributionMapping dmap(ba);
    int ncomps = 6; //  Ex Ey Ez Bx By Bz
    int nghost = 1; //  Ex Ey Ez Bx By Bz
    m_mf_output[0][lev] = amrex::MultiFab( amrex::convert(ba, amrex::IntVect::TheNodeVector()), dmap, ncomps, nghost);
    m_mf_output[0][lev].setVal(0.);
}

void
StationDiagnostics::PrepareFieldDataForOutput ()
{
    const int num_station_functors = 1;
    const int nlev = 1;
    int k_index = m_slice_counter;
    for (int lev = 0; lev < nlev; ++lev) {
        for (int i = 0; i < num_station_functors; ++i) {
            // number of slices = 1
            const bool ZSliceInDomain = GetZSliceInDomain(lev);
            m_all_field_functors[lev][i]->PrepareFunctorData(0, ZSliceInDomain, m_station_loc, m_buffer_box,
                                                             m_slice_counter, m_buffer_size, 0);
        }
    }
}

bool
StationDiagnostics::GetZSliceInDomain (const int lev)
{
    auto & warpx = WarpX::GetInstance();
    const amrex::RealBox& prob_domain = warpx.Geom(lev).ProbDomain();
    if ( ( m_station_loc <= prob_domain.lo(WARPX_ZINDEX) ) or
         ( m_station_loc >= prob_domain.hi(WARPX_ZINDEX) ) )
    {
        return false;
    }
    return true;
}

void
StationDiagnostics::UpdateBufferData ()
{
    if (GetZSliceInDomain(0))
        m_slice_counter++;
    if (m_slice_counter > 0 and !GetZSliceInDomain(0)) {
        m_last_timeslice_filled = true;
    }
}

void
StationDiagnostics::InitializeParticleFunctors ()
{
    auto& warpx = WarpX::GetInstance();
    const int num_station_buffers = 1;
    const MultiParticleContainer& mpc = warpx.GetPartContainer();
    m_all_particle_functors.resize(m_output_species_names.size());
    m_totalParticles_in_buffer[0].resize(m_output_species_names.size());
    for (int i = 0; i < m_all_particle_functors.size(); ++i)
    {
        m_totalParticles_in_buffer[0][i] = 0;
        const int idx = mpc.getSpeciesID(m_output_species_names[i]);
        m_all_particle_functors[i] = std::make_unique<StationParticleFunctor>(mpc.GetParticleContainerPtr(idx), m_output_species_names[i], num_station_buffers);
    }
}

void
StationDiagnostics::InitializeParticleBuffer ()
{
    auto& warpx = WarpX::GetInstance();
    const MultiParticleContainer& mpc = warpx.GetPartContainer();
    const int num_stationdiag_buffers = 1;
    for (int i = 0; i < num_stationdiag_buffers; ++i) {
        m_particles_buffer[i].resize(m_output_species_names.size());
        for (int isp = 0; isp < m_particles_buffer[i].size(); ++isp) {
            m_particles_buffer[i][isp] = std::make_unique<PinnedMemoryParticleContainer>(warpx.GetParGDB());
            const int idx = mpc.getSpeciesID(m_output_species_names[isp]);
            m_output_species[i].push_back(ParticleDiag(m_diag_name,
                                                       m_output_species_names[isp],
                                                       mpc.GetParticleContainerPtr(idx),
                                                       m_particles_buffer[i][isp].get() ));
        }
    }

}

void
StationDiagnostics::PrepareParticleDataForOutput ()
{
    const int lev = 0;
    auto& warpx = WarpX::GetInstance();
    for (int i = 0; i < m_particles_buffer.size(); ++i) {
        for (int isp = 0; isp < m_all_particle_functors.size(); ++isp )
        {
            amrex::Box domain = (warpx.boxArray(lev)).minimalBox();
            domain.setSmall(WARPX_ZINDEX, 0);
            domain.setBig(WARPX_ZINDEX, (m_buffer_size - 1));
            const amrex::Box particle_buffer_box = domain;
            amrex::BoxArray diag_ba;
            diag_ba.define(m_buffer_box);
            amrex::BoxArray ba = diag_ba.maxSize(256);
            amrex::DistributionMapping dmap(ba);
            m_particles_buffer[i][isp]->SetParticleBoxArray(lev, ba);
            m_particles_buffer[i][isp]->SetParticleDistributionMap(lev, dmap);
            m_particles_buffer[i][isp]->SetParticleGeometry(lev, warpx.Geom(0));
        }
    }
    for (int isp = 0; isp < m_all_particle_functors.size(); ++isp)
    {
        m_all_particle_functors[isp]->PrepareFunctorData(0, GetZSliceInDomain(0), m_station_loc);
    }
}

void
StationDiagnostics::Flush (int i_buffer)
{
    if (m_slice_counter == 0) return;
    auto & warpx = WarpX::GetInstance();
    m_tmax = warpx.gett_new(0);
    std::string filename = m_file_prefix;
    filename = amrex::Concatenate(m_file_prefix, i_buffer, m_file_min_digits);
    constexpr int permission_flag_rwxrxrx = 0755;
    if (! amrex::UtilCreateDirectory(filename, permission_flag_rwxrxrx) ) {
        amrex::CreateDirectoryFailed(filename);
    }
    amrex::Print() << " finename : " << filename << "\n";
    WriteStationHeader(filename);
    for (int lev = 0; lev < 1; ++lev) {
        const std::string buffer_path = filename + amrex::Concatenate("/Level_",lev,1) + "/";
        if (m_flush_counter == 0) {
            if (! amrex::UtilCreateDirectory(buffer_path, permission_flag_rwxrxrx) ) {
                amrex::CreateDirectoryFailed(buffer_path);
            }
        }
        std::string buffer_string = amrex::Concatenate("buffer-",m_flush_counter,m_file_min_digits);
        const std::string prefix = amrex::MultiFabFileFullPrefix(lev, filename, "Level_", buffer_string);
        amrex::Print() << " Writing buffer " << buffer_string << "\n";
        amrex::VisMF::Write(m_mf_output[i_buffer][lev], prefix);
    }

    for (int isp = 0; isp < m_particles_buffer[0].size(); ++isp) {
        FlushParticleBuffer(filename, m_output_species_names[isp], isp, i_buffer);
        m_totalParticles_in_buffer[i_buffer][isp] = 0;
    }

    // reset counter
    m_slice_counter = 0;
    // update Flush counter
    m_flush_counter++;
}

void
StationDiagnostics::WriteStationHeader (const std::string& filename)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        const std::string HeaderFileName(filename + "/StationHeader");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);

        HeaderFile.precision(17);

        HeaderFile << m_tmin << "\n";
        HeaderFile << m_tmax << "\n";
        HeaderFile << m_buffer_size << "\n";
        HeaderFile << m_slice_counter << "\n";
        HeaderFile << m_flush_counter << "\n";


    }
}

void
StationDiagnostics::FlushParticleBuffer (std::string path, std::string species_name, int isp, int i_buffer)
{

    amrex::Vector<std::string> real_names;
    amrex::Vector<std::string> int_names;
    amrex::Vector<int> int_flags;
    amrex::Vector<int> real_flags;

    real_names.push_back("weight");
    real_names.push_back("momentum_x");
    real_names.push_back("momentum_y");
    real_names.push_back("momentum_z");

    real_names.resize(m_particles_buffer[i_buffer][isp]->NumRealComps());
    auto runtime_rnames = m_particles_buffer[i_buffer][isp]->getParticleRuntimeComps();
    for (auto const& x : runtime_rnames) {real_names[x.second+PIdx::nattribs] = x.first;}

    real_flags.resize(m_particles_buffer[i_buffer][isp]->NumRealComps(), 1);

    int_names.resize(m_particles_buffer[i_buffer][isp]->NumIntComps());
    auto runtime_inames = m_particles_buffer[i_buffer][isp]->getParticleRuntimeiComps();
    for (auto const& x : runtime_inames) { int_names[x.second + 0] = x.first; }

    int_flags.resize(m_particles_buffer[i_buffer][isp]->NumIntComps(), 1);

    std::string species_dir = path + "/" + species_name;
    // create directory for recorded species
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        if ( ! amrex::UtilCreateDirectory(species_dir, 0755))
        {
            amrex::CreateDirectoryFailed(species_dir);
        }
    }
    amrex::ParallelDescriptor::Barrier();
    // create direcotry for buffer
    std::string buffer_string = amrex::Concatenate("pbuffer-", m_flush_counter, m_file_min_digits);
    std::string buffer_dir = species_dir + "/" + buffer_string;
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        if ( ! amrex::UtilCreateDirectory(buffer_dir, 0755))
        {
            amrex::CreateDirectoryFailed(buffer_dir);
        }
    }
    amrex::ParallelDescriptor::Barrier();

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        WriteHeaderFile(buffer_dir, real_names, int_names, isp);
        WriteParticleData(buffer_dir, real_names, int_names, real_flags, int_flags, isp);
    }
}

void
StationDiagnostics::WriteHeaderFile (std::string pdir, amrex::Vector<std::string> real_names,
                                     amrex::Vector<std::string> int_names, int isp)
{
    std::string HdrFileName = pdir;
    if ( ! HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/') {
        HdrFileName += '/';
    }

    HdrFileName += "Header";

    std::ofstream HdrFile;

    HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);

    if ( ! HdrFile.good()) { amrex::FileOpenFailed(HdrFileName); }

    HdrFile << AMREX_SPACEDIM << "\n";

    // number of real parameters
    HdrFile << real_names.size() << "\n";

    // Real component names
    for (int i = 0; i < real_names.size(); ++i)
    {
        HdrFile << real_names[i] << "\n";
    }

    // number of int parameters
    HdrFile << int_names.size() << "\n";

    // int component names
    for (int i = 0; i < int_names.size(); ++i)
    {
        HdrFile << int_names[i] << "\n";
    }

    // number of particles
    HdrFile << m_totalParticles_in_buffer[0][isp] << "\n";

    // finest level
    HdrFile << m_particles_buffer[0][isp]->finestLevel() << "\n";

    // Box Array size
    for (int lev = 0; lev <= m_particles_buffer[0][isp]->finestLevel(); ++lev)
    {
        HdrFile << m_particles_buffer[0][isp]->ParticleBoxArray(lev).size() << "\n";
    }

    HdrFile.flush();
    HdrFile.close();
    if ( ! HdrFile.good())
    {
        amrex::Abort("ParticleContainer::Checkpoint(): problem writing HdrFile");
    }
}

void
StationDiagnostics::WriteParticleData (std::string pdir, amrex::Vector<std::string> real_names,
                                       amrex::Vector<std::string> int_names,
                                       amrex::Vector<int> real_flags, amrex::Vector<int> int_flags,
                                       int isp)
{
    const int NProcs = amrex::ParallelDescriptor::NProcs();
    const int IOProcNumber = amrex::ParallelDescriptor::IOProcessorNumber();

    // Writing data out in parallel
    // Allow upto nOutFiles activate writers at a time
    int nOutFiles = 256;
    nOutFiles = std::max(1, std::min(nOutFiles, NProcs));

    bool gotsomeParticles = (m_particles_buffer[0][isp] > 0);
    amrex::Print() << " got some particles " << gotsomeParticles << "\n";
    const int lev = 0;

    // flags to determine which particles to output
    amrex::Vector<std::map<std::pair<int,int>, amrex::Vector<int> > > particle_io_flags(m_particles_buffer[0][isp]->GetParticles().size());
    const auto pmap = m_particles_buffer[0][isp]->GetParticles(lev);
    for (const auto & kv : pmap)
    {
        auto& flags = particle_io_flags[lev][kv.first];
        const auto np = kv.second.numParticles(); 
        flags.resize(np,1);
    }
    amrex::Gpu::Device::streamSynchronize();

    amrex::MFInfo info;
    info.SetAlloc(false);
    amrex::MultiFab state(m_particles_buffer[0][isp]->ParticleBoxArray(lev),
                          m_particles_buffer[0][isp]->ParticleDistributionMap(lev),
                          1,0,info);
    // We eventually want to write out the file name and the offset
    // into that file into which each grid of particles is written.
    amrex::Vector<int>  which(state.size(),0);
    amrex::Vector<int > count(state.size(),0);
    amrex::Vector<amrex::Long> where(state.size(),0);

    if (gotsomeParticles) {
        std::string LevelDir = pdir;
        if ( ! LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/') { LevelDir += '/'; }
        LevelDir = amrex::Concatenate(LevelDir.append("Level_"), lev, 1);
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            if ( ! amrex::UtilCreateDirectory(LevelDir, 0755)) {
                amrex::CreateDirectoryFailed(LevelDir);
            }
        }
        amrex::ParallelDescriptor::Barrier();

        std::string filePrefix = LevelDir;
        filePrefix += '/';
        filePrefix += amrex::ParticleContainerBase::DataPrefix();
        bool groupSets(false), setBuf(true);

        for ( amrex::NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf); nfi.ReadyToWrite(); ++nfi)
        {
            auto& myStream = (std::ofstream&) nfi.Stream();
//           pc.WriteParticles(lev, myStream, nfi.FileNumber(), which, count, where,
//                             write_real_comp, write_int_comp, particle_io_flags, is_checkpoint);

            std::map<int, amrex::Vector<int> > tile_map;

            for (const auto& kv : m_particles_buffer[0][isp]->GetParticles(lev))
            {
                const int grid = kv.first.first;
                const int tile = kv.first.second;
                tile_map[grid].push_back(tile);
                const auto flags = particle_io_flags[lev].at(kv.first);

                count[grid] += amrex::particle_detail::countFlags(flags);
            }


            for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi)
            {
                const int grid = mfi.index();

                if (count[grid] == 0) continue;
                amrex::Vector<int> istuff;
                amrex::Vector<amrex::ParticleReal> rstuff;
                amrex::particle_detail::packIOData(istuff, rstuff, m_particles_buffer[0][isp].get(), lev, grid,
                                                   real_flags, int_flags,
                                                   particle_io_flags, tile_map[grid], count[grid], true);

                amrex::writeIntData(istuff.dataPtr(), istuff.size(), myStream);
                myStream.flush();

                amrex::writeDoubleData( (double*) rstuff.dataPtr(), rstuff.size(), myStream);
                myStream.flush();
            }

        }
    }
}
