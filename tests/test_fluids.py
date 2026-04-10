"""Tests for thermosim.fluids — CoolProp property wrappers."""

import pytest

from thermosim.fluids import FluidState, fluid_state, saturation_curve


class TestFluidState:
    def test_fluid_state_from_PT(self):
        """Water at 1 MPa, 450 K — check against known CoolProp values."""
        state = fluid_state("Water", P=1e6, T=450.0)
        assert isinstance(state, FluidState)
        assert state.T == pytest.approx(450.0, rel=1e-6)
        assert state.P == pytest.approx(1e6, rel=1e-6)
        assert state.h == pytest.approx(749.3e3, rel=0.01)
        assert state.s > 0
        assert state.rho > 0

    def test_fluid_state_from_PQ(self):
        """Saturated steam at 100 kPa, quality=1."""
        state = fluid_state("Water", P=100e3, Q=1.0)
        assert state.T == pytest.approx(372.76, rel=0.01)
        assert state.quality == pytest.approx(1.0)

    def test_fluid_state_from_PS(self):
        """Water at 1 MPa, s=6.5 kJ/(kg·K) — superheated region."""
        state = fluid_state("Water", P=1e6, s=6500.0)
        assert state.T > 400
        assert state.P == pytest.approx(1e6, rel=1e-6)

    def test_fluid_state_r134a(self):
        """R-134a at 500 kPa, 300 K."""
        state = fluid_state("R134a", P=500e3, T=300.0)
        assert state.T == pytest.approx(300.0, rel=1e-6)
        assert state.h > 0

    def test_fluid_state_co2(self):
        """CO2 at 5 MPa, 350 K."""
        state = fluid_state("CO2", P=5e6, T=350.0)
        assert state.T == pytest.approx(350.0, rel=1e-6)

    def test_invalid_fluid_raises(self):
        with pytest.raises(ValueError, match="Unknown fluid"):
            fluid_state("FakeFluid123", P=1e6, T=300.0)

    def test_insufficient_inputs_raises(self):
        with pytest.raises(ValueError, match="Exactly two"):
            fluid_state("Water", P=1e6)


class TestSaturationCurve:
    def test_saturation_curve_returns_arrays(self):
        sat = saturation_curve("Water", n_points=50)
        assert len(sat["T"]) == 50
        assert len(sat["s_f"]) == 50
        assert len(sat["s_g"]) == 50
        assert len(sat["h_f"]) == 50
        assert len(sat["h_g"]) == 50
        for sf, sg in zip(sat["s_f"], sat["s_g"]):
            assert sg >= sf

    def test_saturation_curve_temperature_range(self):
        sat = saturation_curve("Water", n_points=100)
        assert sat["T"][0] < 280
        assert sat["T"][-1] > 640
