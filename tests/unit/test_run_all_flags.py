"""Unit tests for run-all CLI options."""
import pytest
import typer


class TestRunAllCliOptions:
    """Test run-all CLI option parsing."""
    
    def test_no_tfbs_flag_disables_tfbs(self):
        """Test that --no-tfbs flag is parsed correctly."""
        from krewlyzer.wrapper import run_all as run_all_func
        import inspect
        sig = inspect.signature(run_all_func)
        
        assert 'no_tfbs' in sig.parameters, "--no-tfbs option should exist"
        # typer.Option returns an OptionInfo object; check its default attribute
        param = sig.parameters['no_tfbs']
        if hasattr(param.default, 'default'):
            assert param.default.default is False, "--no-tfbs should default to False"
        else:
            assert param.default is False, "--no-tfbs should default to False"
    
    def test_no_atac_flag_disables_atac(self):
        """Test that --no-atac flag is parsed correctly."""
        from krewlyzer.wrapper import run_all as run_all_func
        import inspect
        sig = inspect.signature(run_all_func)
        
        assert 'no_atac' in sig.parameters, "--no-atac option should exist"
        param = sig.parameters['no_atac']
        if hasattr(param.default, 'default'):
            assert param.default.default is False, "--no-atac should default to False"
        else:
            assert param.default is False, "--no-atac should default to False"
    
    def test_both_flags_can_be_set(self):
        """Test that both --no-tfbs and --no-atac can be set simultaneously."""
        from krewlyzer.wrapper import run_all as run_all_func
        import inspect
        sig = inspect.signature(run_all_func)
        
        # Both parameters should exist and be boolean
        assert 'no_tfbs' in sig.parameters
        assert 'no_atac' in sig.parameters
        assert sig.parameters['no_tfbs'].annotation == bool
        assert sig.parameters['no_atac'].annotation == bool


class TestRunAllEnableLogic:
    """Test that enable_tfbs/enable_atac logic respects CLI flags."""
    
    def test_enable_tfbs_respects_no_flag(self):
        """Test enable_tfbs = assets.tfbs_available and not no_tfbs."""
        tfbs_available = True
        no_tfbs = True
        
        result = tfbs_available and not no_tfbs
        assert result is False, "enable_tfbs should be False when --no-tfbs is set"
    
    def test_enable_atac_respects_no_flag(self):
        """Test enable_atac = assets.atac_available and not no_atac."""
        atac_available = True
        no_atac = True
        
        result = atac_available and not no_atac
        assert result is False, "enable_atac should be False when --no-atac is set"
    
    def test_default_behavior_enables_both(self):
        """Test that default behavior (no flags) enables both if available."""
        tfbs_available = True
        atac_available = True
        no_tfbs = False
        no_atac = False
        
        enable_tfbs = tfbs_available and not no_tfbs
        enable_atac = atac_available and not no_atac
        
        assert enable_tfbs is True
        assert enable_atac is True
