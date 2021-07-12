"""
htc_fs:
TODO: we will use a consant heat transfer coefficient here, which is 
not ideal for the future, because the heat transfer coefficient depends 
on Re and Pr and therefore is dependet of flow speed (turbulence) and 
the thermal properties near the "wall".
Therefore, we will have to expand this function to be dependend of 
(Tf,Tg,material_solid, material_fluid) in the future ..
"""
function htc_fs(ht_type::ConstantHeatTransfer)
    return ht_type.val
end