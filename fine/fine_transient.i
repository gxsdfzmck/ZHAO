# Flow and solute transport along a fracture embedded in a porous matrix
# The fracture is represented by lower dimensional elements
# fracture aperture = 6e-4m
# fracture porosity = 6e-4m = phi * a
# fracture permeability = 1.8e-11 which is based on k=3e-8 from a**2/12, and k*a = 3e-8*6e-4
# matrix porosity = 0.1
# matrix permeanility = 1e-20

[Mesh]
  file = 'steady_out_fine.e' #'steady_out_fine.e'
  block_id = '1 2 3'
  block_name = 'fracture matrix1 matrix2'
 #
 # boundary_id = '1 2 1111 2222 3333 4444 5555 6666'
 # boundary_name = 'inject produce top bottom left right front back'
[]

[GlobalParams]
  displacements = 'sdisp_x sdisp_y sdisp_z'
  PorousFlowDictator = dictator
  gravity = '0 0 -9.8'
[]


[Variables]
  [./pp]
    scaling = 1E9
  [../]
  [./T]
    scaling = 1E3
  [../]
  [./sdisp_x]
  [../]
  [./sdisp_y]
  [../]
  [./sdisp_z]
  [../]
[]

[AuxVariables]
  [./velocity_x]
    family = MONOMIAL
    order = CONSTANT
    block = 'fracture'
  [../]
  [./velocity_y]
    family = MONOMIAL
    order = CONSTANT
    block = 'fracture'
  [../]
  [./velocity_z]
    family = MONOMIAL
    order = CONSTANT
    block = 'fracture'
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./velocity_x]
    type = PorousFlowDarcyVelocityComponentLowerDimensional
    variable = velocity_x
    component = x
    aperture = 5E-5            
  [../]
  [./velocity_y]
    type = PorousFlowDarcyVelocityComponentLowerDimensional
    variable = velocity_y
    component = y
    aperture = 5E-5
  [../]
  [./velocity_z]
    type = PorousFlowDarcyVelocityComponentLowerDimensional
    variable = velocity_z
    component = z
    aperture = 5E-5
  [../]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./stress_xz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 0
    index_j = 2
  [../]
  [./stress_yx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yx
    index_i = 1
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
  [../]
  [./stress_zx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zx
    index_i = 2
    index_j = 0
  [../]
  [./stress_zy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zy
    index_i = 2
    index_j = 1
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
[]

[ICs]
  [./pp_matrix]
   type = FunctionIC 
   function = initial_pp
   variable = pp
  [../]
  [./T_matrix]
   type = FunctionIC 
   function = initial_T
   variable = T
  [../]
[]

[Functions]
  [./initial_pp]
    type = SolutionFunction
    from_variable = pp
    solution = steady_solution_pp
  [../]
  [./initial_T]
    type = SolutionFunction
    from_variable = T
    solution = steady_solution_T
  [../]
  [./top_force_pressure]
    type = ParsedFunction
    value = '1E8'
  [../]
  [./front_force_pressure]
    type = ParsedFunction
    value = '80+(5000-z)*(180-80)/5000'
  [../]
  [./right_force_pressure]
    type = ParsedFunction
    value = '120+(5000-z)*(270-120)/5000'
  [../]
[]

[BCs]
  [./ptop]
    type = NeumannBC 
    variable = pp
    boundary =  top
    value = 0
  [../]
  [./pbottom]
    type = NeumannBC 
    variable = pp
    boundary = bottom 
    value = 0 
  [../]
  [./pleft]
    type = NeumannBC 
    variable = pp
    boundary = left
    value = 0
  [../]
  [./pright]
    type = NeumannBC 
    variable = pp
    boundary = right 
    value = 0
  [../]
  [./pfront]
    type = NeumannBC 
    variable = pp
    boundary = front
    value = 0
  [../]
  [./pback]
    type = NeumannBC 
    variable = pp
    boundary = back
    value = 0
  [../]
  [./pInject]
    type = PresetBC
    variable = pp
    boundary = inject
    value = 27.3E6
  [../]
  [./pProduce]
    type = PresetBC
    variable = pp
    boundary = produce
    value = 9.8E6
  [../]
  [./Tinject]
    type = PresetBC
    variable = T
    boundary = inject
    value = 303    
  [../]
  ##### Energy BC #######
  [./Ttop]
    type = NeumannBC
    variable = T
    boundary = top
    value = 0
  [../]
  [./Tbottom]
    type = NeumannBC
    variable = T
    boundary = bottom
    value = 0
  [../]
  [./Tleft]
    type = NeumannBC
    variable = T
    boundary = left
    value = 0
  [../]
  [./Tright]
    type = NeumannBC
    variable = T
    boundary = right
    value = 0
  [../]
  [./Tfront]
    type = NeumannBC
    variable = T
    boundary = front
    value = 0
  [../]
  [./Tback]
    type = NeumannBC
    variable = T
    boundary = back
    value = 0
  [../]
  [./solid_top]
    type = Pressure
    variable = sdisp_z
    component = 2
    function = top_force_pressure  # 100MPa
    boundary = top
  [../]
########## solid mechanics BCS ###########
  [./solid_bottom]
    type = PresetBC
    variable = sdisp_z
    value = 0
    boundary = bottom
  [../]
  [./solid_left]
    type = PresetBC
    variable = sdisp_y
    value = 0
    boundary = left
  [../]
  [./solid_right]
    type = Pressure
    variable = sdisp_y
    component = 1
    function = right_force_pressure
    use_displaced_mesh = false
    boundary = right
  [../]
  [./solid_front]
    type = Pressure
    variable = sdisp_x
    component = 0
    use_displaced_mesh = false
    function = front_force_pressure
    boundary = front
  [../]
  [./solid_back]
    type = PresetBC
    variable = sdisp_x
    value =  0
    boundary = back
  [../]
[]

[Kernels]
  [./massDeriv]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0 
    variable = pp
  [../]
  [./adv0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pp
  [../]
 ############ Energy Kernels for matrix ##########
  [./EnergyTimeDeriv_matrix]
    type = PorousFlowEnergyTimeDerivative
    variable = T
  [../]
  [./EnergyAdvection_matrix]
    type = PorousFlowHeatAdvection
    fluid_component = 0
    variable = T
    gravity = '0 0 -9.8'
  [../]
  [./EnergyConduciton]
    type = PorousFlowHeatConduction
    variable = T
  [../]
######### Tensor Mechanicals Kernels ########
  [./grad_stress_x]
    type = StressDivergenceTensors
    variable = sdisp_x
    component = 0
  [../]
  [./grad_stress_y]
    type = StressDivergenceTensors
    variable = sdisp_y
    component = 1
  [../]
  [./grad_stress_z]
    type = StressDivergenceTensors
    variable = sdisp_z
    component = 2
  [../]
  [./gravity_z]
    type = Gravity
    variable = sdisp_z
    value = -9.8
  [../]
  [./poro_x]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 1.0
    variable = sdisp_x
    component = 0
  [../]
  [./poro_y]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 1.0
    variable = sdisp_y
    component = 1
  [../]
  [./poro_z]
    type = PorousFlowEffectiveStressCoupling
    biot_coefficient = 1.0
    variable = sdisp_z
    component = 2
  [../]
  [./poro_vol_exp]
    type = PorousFlowMassVolumetricExpansion
    variable = pp
    fluid_component = 0
  [../]
  [./heat_vol_exp]
    type = PorousFlowHeatVolumetricExpansion
    variable = T
  [../]
[]

[UserObjects]
  [./steady_solution_pp]
    type = SolutionUserObject
    timestep = LATEST
    system_variables = 'pp' 
    mesh = steady_out_fine.e
  [../]
  [./steady_solution_T]
    type = SolutionUserObject
    timestep = LATEST
    system_variables = 'T' 
    mesh = steady_out_fine.e
  [../]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp T sdisp_x sdisp_y sdisp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  [../]
[]

[Modules]
  [./FluidProperties]
   # [./simple_fluid]
   #   type = SimpleFluidProperties
   #   bulk_modulus = 2e9
   #   density0 = 1000
   #   thermal_expansion = 0
   #   viscosity = 1e-3
   # [../]
    [./water97property]
      type = Water97FluidProperties
    [../]
  [../]
[]

[Materials]
  [./temperature_matrix]
    type = PorousFlowTemperature
    temperature = T
  [../]
  [./temperature_nodal_matrix]
    type = PorousFlowTemperature
    at_nodes = true
    temperature = T
  [../]
  [./rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 1.08E3
    density = 2.7E3
  [../]
  [./ppss]
    type = PorousFlow1PhaseFullySaturated
    at_nodes = true
    porepressure = pp
  [../]
  [./ppss_qp]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  [../]
  [./water97property]
    type = PorousFlowSingleComponentFluid
    fp = water97property
    phase = 0
    at_nodes = true
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
    at_nodes = true
  [../]
  [./simple_fluid_qp]
    type = PorousFlowSingleComponentFluid
    fp = water97property 
    phase = 0
  [../]
  [./thermal_conductivity_matrix]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '3.0 0 0 0 3.0 0 0 0 3.0'
    block = 'matrix1 matrix2'
  [../]
  [./thermal_conductivity_fracture]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '1.0 0 0 0 1.0 0 0 0 1.0'
    block = 'fracture'
  [../]
  [./poro_fracture]
    type = PorousFlowPorosity
    PorousFlowDictator = dictator
    thermal = true
    ensure_positive = true
    fluid = true
    mechanical = true
    biot_coefficient = 1.0
    porosity_zero = 0.1
    solid_bulk = 2.30E10
    thermal_expansion_coeff = 1.05E-5
   # reference_porepressure = 'pp'
   # reference_temperature = 'T'
    block = 'fracture'
  [../]
  [./poro_matrix]
    type = PorousFlowPorosity
    PorousFlowDictator = dictator
    ensure_positive = true
    thermal = true
    fluid = true
    mechanical = true
    porosity_zero = 0.01
   # reference_porepressure = 'pp'
   # reference_temperature = 'T'
    solid_bulk = 2.30E10 
    biot_coefficient = 1
    thermal_expansion_coeff = 1.0E-5
    block = 'matrix1 matrix2'
  [../]
  [./poro_fracture_nodal]
    at_nodes = true
    type = PorousFlowPorosity
    PorousFlowDictator = dictator
    ensure_positive = true
    thermal = true
    fluid = true
    mechanical = true
    biot_coefficient = 1.0
    porosity_zero = 0.1
    solid_bulk = 2.30E10
    thermal_expansion_coeff = 1.05E-5
   # reference_porepressure = 'pp'
   # reference_temperature = 'T'
    block = 'fracture'
  [../]
  [./poro_matrix_nodal]
    at_nodes = true
    type = PorousFlowPorosity
    PorousFlowDictator = dictator
    ensure_positive = true
    thermal = true
    fluid = true
    mechanical = true
    porosity_zero = 0.01
   # reference_porepressure = 'pp'
   # reference_temperature = 'T'
    solid_bulk = 2.30E10 
    biot_coefficient = 1
    thermal_expansion_coeff = 1.0E-5
    block = 'matrix1 matrix2'
  [../]
  [./diff1]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '1e-9'
    tortuosity = 1.0
    block = 'fracture'
  [../]
  [./diff2]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '1e-9'
    tortuosity = 0.1
    block = 'matrix1 matrix2'
  [../]
  [./permeability_fracture]
   # type = PorousFlowPermeabilityConst
   # permeability = '1.05e-8 0 0 0 1.05e-8 0 0 0 1.05e-8'   # 1.8e-11 = a * kf
    type = PorousFlowPermeabilityKozenyCarman
    poroperm_function = kozeny_carman_phi0
    PorousFlowDictator = dictator
    m = 2
    phi0 = 0.25
    k0 = 1.1E-8
    n = 3
    block = 'fracture'
  [../]
  [./permeability_matrix]
   # type = PorousFlowPermeabilityConst
   # permeability = '1.096e-16 0 0 0 1.096e-16 0 0 0 1.096e-16'
    type = PorousFlowPermeabilityKozenyCarman
    poroperm_function = kozeny_carman_phi0
    PorousFlowDictator = dictator
    m = 2
    phi0 = 0.15
    k0 = 1.096E-16
    n = 3
    block = 'matrix1 matrix2'
  [../]
  [./relp]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  [../]
  [./relp_nodal]
    type = PorousFlowRelativePermeabilityConst
    at_nodes = true
    phase = 0
  [../]
########## solid mechanics material########
  [./elastic_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 4.84E10
    poissons_ratio = 0.15
   # type = ComputeElasticityTensor
   # C_ijkl = '3.478E10 2.10E9' # young = 48.4GPa, poisson = 0.15
   # fill_method = symmetric_isotropic
  [../]
  [./thermal_expansion_strain]
    type = ComputeThermalExpansionEigenstrain
    stress_free_temperature = 293
    thermal_expansion_coeff = 1.0E-5
    temperature = T
    eigenstrain_name = eigenstrain
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'sdisp_x sdisp_y sdisp_z'
   # type = ComputeIncrementalSmallStrain
  [../]
  [./density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2.683E3 # (1-0.01)*2700 + 0.01*1.0E3
  [../]
  [./eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    at_nodes = true
  [../]
  [./eff_fluid_pressure_qp]
    type = PorousFlowEffectiveFluidPressure
  [../]
  [./stress]
    type = ComputeLinearElasticStress
   # type = ComputeMultiPlasticityStress
   # ep_plastic_tolerance = 1E-5
  [../]
  [./vol_strain]
    type = PorousFlowVolumetricStrain
  [../]
[]

[Preconditioning]
  [./basic]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2             '
  [../]
[]

[Executioner]
 # type = Transient
 # num_steps = 20
 # solve_type = NEWTON
 # dt = 86400.0

 # l_max_its = 20
 # l_tol = 1e-5
 # nl_max_its = 30
 # nl_abs_tol = 1e-9
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_hypre_type

                         -ksp_gmres_restart -snes_ls

                         -pc_hypre_boomeramg_strong_threshold'

  petsc_options_value = 'hypre boomeramg 201 cubic 0.7'
  
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 8640
  [../]
    end_time = 4.32E7 # 500 day 
    l_tol = 1e-7
    l_max_its = 500
    nl_rel_tol = 1e-10
    nl_abs_tol = 1e-7
[]

#[VectorPostprocessors]
#  [./x_pp]
#    type = LineValueSampler
#    start_point = '0 0 0'
#    end_point = '1000 0 0'
#    sort_by = x
#    num_points = 50
#    variable = pp
#    outputs = csv
#  [../]
#[]

[Outputs]
 # [./csv]
 #   type = CSV
 #   execute_on = 'final'
 # [../]
  [./exduos]
    type = Exodus
    file_base = transient_out_fine
    execute_on = 'timestep_end'
  [../]
[]
