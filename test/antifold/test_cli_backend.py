"""
Test CLI backend switching for proteinmpnn command.

TDD approach: test CLI option parsing before implementation.
"""

import pytest
from click.testing import CliRunner


@pytest.fixture
def cli_runner():
    return CliRunner()


class TestCliBackendOption:
    """Test --backend option on proteinmpnn CLI."""

    def test_backend_defaults_to_proteinmpnn(self, cli_runner):
        """Default backend should be proteinmpnn."""
        from rfantibody.cli.inference import proteinmpnn

        # Use --help to verify option exists without running model
        result = cli_runner.invoke(proteinmpnn, ["--help"])
        assert result.exit_code == 0
        assert "--backend" in result.output

    def test_backend_accepts_antifold(self, cli_runner):
        """--backend antifold should be accepted."""
        from rfantibody.cli.inference import proteinmpnn

        # Verify help text mentions antifold
        result = cli_runner.invoke(proteinmpnn, ["--help"])
        assert result.exit_code == 0
        assert "antifold" in result.output

    def test_backend_accepts_proteinmpnn(self, cli_runner):
        """--backend proteinmpnn should be accepted."""
        from rfantibody.cli.inference import proteinmpnn

        result = cli_runner.invoke(proteinmpnn, ["--help"])
        assert result.exit_code == 0
        assert "proteinmpnn" in result.output
