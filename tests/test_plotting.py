"""Tests for thermosim.plotting — diagram generators."""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest
from thermosim.plotting import plot_ts_diagram, plot_ph_diagram, apply_style


class TestApplyStyle:
    def test_applies_without_error(self):
        fig, ax = plt.subplots()
        apply_style(ax, title="Test", xlabel="X", ylabel="Y")
        assert ax.get_title() == "Test"
        assert ax.get_xlabel() == "X"
        plt.close(fig)


class TestTSDiagram:
    def test_returns_figure_and_axes(self):
        fig, ax = plot_ts_diagram("Water")
        assert isinstance(fig, plt.Figure)
        assert isinstance(ax, plt.Axes)
        assert len(ax.lines) >= 2
        plt.close(fig)

    def test_with_state_points(self):
        from thermosim.fluids import fluid_state
        states = [
            fluid_state("Water", P=10e3, Q=0.0),
            fluid_state("Water", P=10e3, Q=1.0),
        ]
        fig, ax = plot_ts_diagram("Water", states=states)
        assert len(ax.lines) >= 2
        plt.close(fig)


class TestPHDiagram:
    def test_returns_figure_and_axes(self):
        fig, ax = plot_ph_diagram("Water")
        assert isinstance(fig, plt.Figure)
        assert isinstance(ax, plt.Axes)
        assert len(ax.lines) >= 2
        plt.close(fig)
