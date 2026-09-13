/* Copyright 2026 The WarpX Community
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_MLMG_TUNING_CACHE_H_
#define ABLASTR_MLMG_TUNING_CACHE_H_

#include <AMReX.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_FabFactory.H>
#include <AMReX_FArrayBox.H>

#include <map>
#include <memory>
#include <utility>

namespace ablastr::fields
{

/** Tuned MLMG parameters, cached per (level, BoxArray, DistributionMapping).
 *
 * The key holds raw reference ids, so an entry must be erased before its
 * BoxArray or DistributionMapping is replaced (e.g. on load balance).
 */
struct MLMGTuningCache
{
    struct Parms {
        bool single_box = false;
        int mcl = 30; // default in amrex
        int ags = -1; // default in amrex
        bool no_gpu_sync = true; // GPU only: single stream and no stream sync
    };

    using Key = std::pair<int, amrex::FabArrayBase::BDKey>;

    //! Enable/disable tuning (opt-in, since it performs trial solves).
    static void setTuning (bool enable) { m_tune = enable; }
    static bool tuning () { return m_tune; }

    static Key makeKey (int lev, amrex::BoxArray const& ba,
                        amrex::DistributionMapping const& dm)
    {
        return {lev, amrex::FabArrayBase::BDKey(ba.getRefID(), dm.getRefID())};
    }

    //! Returns nullptr if there is no entry.
    static Parms const* find (Key const& key)
    {
        auto const it = m_parms_map.find(key);
        return (it != m_parms_map.end()) ? &(it->second) : nullptr;
    }

    static void insert (Key const& key, Parms const& parms)
    {
        if (m_parms_map.size() > 1000) { m_parms_map.clear(); } // just in case
        m_parms_map[key] = parms;
    }

    //! Erase the entry of a level's current layout before it is replaced.
    static void erase (int lev, amrex::BoxArray const& ba,
                       amrex::DistributionMapping const& dm)
    {
        m_parms_map.erase(makeKey(lev, ba, dm));
    }

    //! Layout of the single-Box solve, plus its EB factory (built by the caller on demand).
    struct SingleBox {
        amrex::Box domain;
        amrex::BoxArray ba;
        amrex::DistributionMapping dm;
        std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>> eb_factory;
    };

    //! Single-Box layout of a level. Rebuilt if the domain or index type changed.
    static SingleBox& singleBox (int lev, amrex::Box const& domain, amrex::IndexType typ)
    {
        if (m_single_box_map.empty()) {
            amrex::ExecOnFinalize([] () { m_single_box_map.clear(); });
        }
        auto& sb = m_single_box_map[lev];
        if (sb.ba.empty() || sb.domain != domain || sb.ba.ixType() != typ) {
            sb = SingleBox{};
            sb.domain = domain;
            sb.ba = amrex::BoxArray(domain);
            sb.ba.convert(typ);
            sb.dm = amrex::DistributionMapping(amrex::Vector<int>{0});
        }
        return sb;
    }

private:
    inline static bool m_tune = false;
    inline static std::map<Key, Parms> m_parms_map;
    inline static std::map<int, SingleBox> m_single_box_map;
};

} // namespace ablastr::fields

#endif // ABLASTR_MLMG_TUNING_CACHE_H_
