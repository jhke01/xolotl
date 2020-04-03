#pragma once

#include <experimental/ReactionNetwork.h>
#include <experimental/NEReaction.h>
#include <experimental/NETraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
class NEReactionGenerator;
}

class NEReactionNetwork : public ReactionNetwork<NEReactionNetwork>
{
    friend class ReactionNetwork<NEReactionNetwork>;
    friend class detail::ReactionNetworkWorker<NEReactionNetwork>;

public:
    using Superclass = ReactionNetwork<NEReactionNetwork>;
    using Subpaving = typename Superclass::Subpaving;
    using Composition = typename Superclass::Composition;
    using Species = typename Superclass::Species;
    using IndexType = typename Superclass::IndexType;
    using ConcentrationsView = typename Superclass::ConcentrationsView;
    using FluxesView = typename Superclass::FluxesView;

    using Superclass::Superclass;

private:
    double
    checkLatticeParameter(double latticeParameter)
    {
        if (latticeParameter <= 0.0) {
            return uraniumDioxydeLatticeConstant;
        }
        return latticeParameter;
    }

    double checkImpurityRadius(double impurityRadius)
    {
        if (impurityRadius <= 0.0) {
            return xenonRadius;
        }
        return impurityRadius;
    }

    void
    addReactionFluxes(ConcentrationsView concentrations, FluxesView fluxes,
        IndexType gridIndex)
    {
        // Does nothing
    }

    void
    addReactionPartials(ConcentrationsView concentrations,
        Kokkos::View<double*> values, IndexType gridIndex)
    {
        // Does nothing
    }

    detail::NEReactionGenerator
    getReactionGenerator() const noexcept;
};

namespace detail
{
class NEReactionGenerator :
    public ReactionGenerator<NEReactionNetwork, NEReactionGenerator>
{
public:
    using Network = NEReactionNetwork;
    using Subpaving = typename Network::Subpaving;
    using Superclass =
        ReactionGenerator<NEReactionNetwork, NEReactionGenerator>;

    using Superclass::Superclass;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    operator()(IndexType i, IndexType j, TTag tag) const;
};
}
}
}

#include <experimental/NEClusterGenerator.h>
#include <experimental/NEReactionNetwork.inl>
