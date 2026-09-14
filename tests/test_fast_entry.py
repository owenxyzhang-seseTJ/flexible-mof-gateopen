import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from run_fast_screen import fast_spec


def base():
    return dict(phases={"closed":{},"open":{}},guest=dict(name="CO2",smiles="O=C=O"))


def test_fast_defaults_and_no_slow_route(tmp_path):
    spec=base();resolved=fast_spec(spec,tmp_path)
    assert resolved["host"]["method"]=="elastic-network"
    assert resolved["host"]["calculator"]["model"]==str(tmp_path/"assets/models/mofs_v2.model")
    assert resolved["raspa"]["seeds"]==[417,941]
    assert "host" not in spec
    for host in [dict(method="phonons"),dict(relax=True),dict(calculator=dict(backend="chgnet"))]:
        with pytest.raises(ValueError):fast_spec(dict(base(),host=host),tmp_path)
    spec=base();spec["guest"]["refinement"]="gfn2-xtb"
    with pytest.raises(ValueError,match="electronic"):fast_spec(spec,tmp_path)


def test_fast_zero_branch_requires_sensitivity(tmp_path):
    spec=base();spec["adsorption_models"]={"closed":{"kind":"nonadsorbing"}}
    with pytest.raises(ValueError,match="sensitivity"):fast_spec(spec,tmp_path)
    spec["adsorption_models"]["closed"]["sensitivity_K_molecules_fu_bar"]=[1e-6]
    assert fast_spec(spec,tmp_path)["adsorption_models"]==spec["adsorption_models"]
