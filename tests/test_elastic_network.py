import numpy as np
import pytest
from ase import Atoms
from flexible_mof_gateopen.elastic_network import network_thermodynamics,compare_states


def test_equivalent_supercell_entropy_and_free_energy():
    # Synthetic periodic spring crystal: a numerical test, not a MOF model.
    a=Atoms("C",positions=[[0,0,0]],cell=[1.5,1.5,1.5],pbc=True)
    b=a.repeat((2,1,1))
    p={"packing_k":0}
    x=network_thermodynamics(a,1,268,mesh=(4,4,4),parameters=p)
    y=network_thermodynamics(b,2,268,mesh=(2,4,4),parameters=p)
    for key in ("S_proxy_J_mol_fu_K","U_vib_proxy_kJ_mol_fu","F_vib_proxy_kJ_mol_fu"):
        assert x[key]==pytest.approx(y[key],abs=1e-8)
    assert x["F_vib_proxy_kJ_mol_fu"]==pytest.approx(x["U_vib_proxy_kJ_mol_fu"]-268*x["S_proxy_J_mol_fu_K"]/1000)


def test_rotation_and_permutation_invariance():
    a=Atoms("C2",positions=[[0,0,0],[1.5,0,0]],cell=[3,1.5,1.5],pbc=True)
    b=a[[1,0]];b.rotate(37,"z",rotate_cell=True)
    r=compare_states([a,b],[0,1],[2,2],268,mesh=(2,2,2),parameter_sets=[{"name":"test","packing_k":0}])
    assert r["parameter_scenarios"][0]["delta_S_proxy_J_mol_fu_K"]==pytest.approx(0,abs=1e-8)
    assert r["parameter_scenarios"][0]["delta_F_host_proxy_kJ_mol_fu"]==pytest.approx(1,abs=1e-8)


def test_proxy_rejects_bad_connectivity_parameters_and_composition():
    a=Atoms("C",positions=[[0,0,0]],cell=[10,10,10],pbc=True)
    with pytest.raises(ValueError,match="unconnected"):
        network_thermodynamics(a,1,268)
    with pytest.raises(ValueError,match="unknown"):
        network_thermodynamics(a,1,268,parameters={"made_up":1})
    b=Atoms("N",positions=[[0,0,0]],cell=[1.5]*3,pbc=True)
    with pytest.raises(ValueError,match="compositions"):
        compare_states([a,b],[0,0],[1,1],268)
