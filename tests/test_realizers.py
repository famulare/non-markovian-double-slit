import numpy as np

from non_markovian_double_slit.finite import FiniteDoubleSlitModel
from non_markovian_double_slit.realizers import (
    IndependentPathRealizer,
    TransportPathRealizer,
    empirical_marginals,
    monotone_transport_coupling,
)


def test_monotone_coupling_has_requested_marginals() -> None:
    source = np.array([0.2, 0.5, 0.3])
    target = np.array([0.4, 0.1, 0.5])
    coupling = monotone_transport_coupling(source, target)
    assert np.all(coupling >= 0)
    assert np.allclose(coupling.sum(axis=0), source)
    assert np.allclose(coupling.sum(axis=1), target)


def test_two_realizers_share_marginals_but_not_path_roughness() -> None:
    model = FiniteDoubleSlitModel.create(
        size=31, length=31.0, slit_separation=6.0, slit_width=2.0
    )
    times = (0.0, 0.7, 1.4)
    process = model.process(times)
    independent = IndependentPathRealizer(process, times)
    transport = TransportPathRealizer(process, times)
    transport_cost = transport.expected_squared_step(model.grid)
    independent_cost = independent.expected_squared_step(model.grid)
    assert transport_cost < independent_cost

    expected = np.vstack(
        [process.initial_distribution, *(process.marginal(time) for time in times)]
    )
    for seed, realizer in enumerate((independent, transport)):
        paths = realizer.sample(30_000, rng=np.random.default_rng(seed))
        observed = empirical_marginals(paths, process.size)
        assert np.max(np.abs(observed - expected)) < 0.015


def test_upper_given_right_can_span_the_full_unit_interval() -> None:
    model = FiniteDoubleSlitModel.create(
        size=51, length=51, slit_separation=8, slit_width=3
    )
    p, r = np.abs(model.slit_state) ** 2, model.probability(4)
    upper, right = model.upper_mode != 0, model.grid > 0
    monotone = monotone_transport_coupling(p, r)
    reversed_rank = monotone_transport_coupling(p[::-1], r)[:, ::-1]
    for joint, expected in ((monotone, 1), (reversed_rank, 0), (np.outer(r, p), 0.5)):
        np.testing.assert_allclose(joint.sum(axis=0), p, atol=1e-14)
        np.testing.assert_allclose(joint.sum(axis=1), r, atol=1e-14)
        assert np.isclose(joint[np.ix_(right, upper)].sum() / r[right].sum(), expected)
