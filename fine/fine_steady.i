# Flow and solute transport along a fracture embedded in a porous matrix
# The fracture is represented by lower dimensional elements
# fracture aperture = 6e-4m
# fracture porosity = 6e-4m = phi * a
# fracture permeability = 1.8e-11 which is based on k=3e-8 from a**2/12, and k*a = 3e-8*6e-4
# matrix porosity = 0.1
# matrix permeanility = 1e-20

[Mesh]
  type = FileMesh
  file = 'mesh_fine.e'
  block_id = '1 2 3'
  block_name = 'fracture matrix1 matrix2'

  boundary_id = '1 2 1111 2222 3333 4444 5555 6666'
  boundary_name = 'inject produce top bottom left right front back'
[]

[GlobalParams]
  displacements = 'sdisp_x sdisp_y sdisp_z'
  PorousFlowDictator = dictator
[]

[Variables]
  [./pp]
    scaling = 1.0E-14
  [../]

  [./T]
    scaling = 1.0E-14
  [../]

  [./sdisp_x]
    scaling = 1.0E-15
  [../]

  [./sdisp_y]
    scaling = 1.0E-15
  [../]

  [./sdisp_z]
    scaling = 1.0E-15
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
    gravity = '0 0 -9.8'
  [../]

  [./velocity_y]
    type = PorousFlowDarcyVelocityComponentLowerDimensional
    variable = velocity_y
    component = y
    aperture = 5E-5
    gravity = '0 0 -9.8'
  [../]

  [./velocity_z]
    type = PorousFlowDarcyVelocityComponentLowerDimensional
    variable = velocity_z
    component = z
    aperture = 5E-5
    gravity = '0 0 -9.8'
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
    variable = pp
    function = '1.013E5+1.0E3*9.8*4000-1.0E3*9.8*z'
  [../]

  [./T_initial]
    type = FunctionIC
    variable = T
    function = '298.15+4000*50/1000-50*z/1000'    
  [../]  
[]

[Functions]
  [./weight_fcn]
    type = ParsedFunction
    vars = 'g p0 rho0 biot'
    vals = '9.8 101325 1000 1.0'
    value = 'biot*(p0-rho0*g*z)-(1.0E8-2683*g*z)'
  [../]   
 
  [./kxx_fcn]
    type = ParsedFunction
    vars = 'g B p0 rho0 biot'
    vals = '9.8 2.0E9 101325 1000 1.0'
    value = 'biot*(101325-rho0*g*z)-(80*1.0E6-(180-80)/5000*1.0E6*z)'
  [../]
  
  [./kyy_fcn]
    type = ParsedFunction
    vars = 'g B p0 rho0 biot'
    vals = '9.8 2.0E9 101325 1000 1.0'
    value = 'biot*(101325-rho0*g*z)-(120*1.0E6-(270-120)/5000*1.0E6*z)'
  [../]
  
[]


[BCs]
  ### pore pressure BC ###
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
  ##### Energy BC #######
  [./Ttop]
    type = FunctionPresetBC 
    variable = T
    boundary = top
    function = '298.15+4000*50/1000-50*z/1000'    
  [../]

  [./Tbottom]
    type = FunctionPresetBC
    variable = T
    boundary = bottom
    function = '298.15+4000*50/1000-50*z/1000'    
  [../]

  [./Tleft]
    type = FunctionPresetBC 
    variable = T
    boundary = left
    function = '298.15+4000*50/1000-50*z/1000'    
  [../]

  [./Tright]
    type = FunctionPresetBC 
    variable = T
    boundary = right
    function = '298.15+4000*50/1000-50*z/1000'    
  [../]

  [./Tfront]
    type = FunctionPresetBC 
    variable = T
    boundary = front
    function = '298.15+4000*50/1000-50*z/1000'    
  [../]

  [./Tback]
    type = FunctionPresetBC 
    variable = T
    boundary = back
    function = '298.15+4000*50/1000-50*z/1000'    
  [../]
 # [./pInject]
 #   type = PresetBC
 #   variable = pp
 #   boundary = inject
 #   value = 27.3E6
 # [../]
 # [./pProduce]
 #   type = PresetBC
 #   variable = pp
 #   boundary = produce
 #   value = 9.8E6
 # [../]
  [./solid_top]
    type = Pressure
    variable = sdisp_z
    component = 2
    factor = 1.0E8
    boundary = top
  [../]

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
    function = '(120-z*(270-120)/5000)*1.0E6'  
    use_displaced_mesh = false
    boundary = right
  [../]

  [./solid_front]
    type = Pressure
    variable = sdisp_x
    component = 0
    function = '(80-z*(180-80)/5000)*1.0E6'
    use_displaced_mesh = false
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
####### Pore press Kernels ########
 # [./massDeriv]
 #   type = PorousFlowMassTimeDerivative
 #   fluid_component = 0
 #   variable = pp
 # [../]
  [./adv0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pp
    gravity = '0 0 -9.8'
  [../]
 ############ Energy Kernels ##########
 # [./EnergyAdvection]
 #   type = PorousFlowHeatAdvection
 #   fluid_component = 0
 #   variable = T
 #   gravity = '0 0 -9.8'
 # [../]
  [./EnergyConduciton]
    type = PorousFlowHeatConduction
    variable = T
  [../]
##### solid mechancis #########
  [./grad_stress_x]
    type = StressDivergenceTensors
    variable = sdisp_x
    component = 0
    use_displaced_mesh = false
  [../]
  [./grad_stress_y]
    type = StressDivergenceTensors
    variable = sdisp_y
    component = 1
    use_displaced_mesh = false
  [../]
  [./grad_stress_z]
    type = StressDivergenceTensors
    variable = sdisp_z
    component = 2
    use_displaced_mesh = false
  [../]
  [./gravity_z]
    type = Gravity
    variable = sdisp_z
    value = -9.8
    use_displaced_mesh = false
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
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp T sdisp_x sdisp_y sdisp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  [../]
  [./ts]
    # large so that there is no plastic deformation
    type = TensorMechanicsHardeningConstant
    value = 1E16
  [../]
[]

[Modules]
  [./FluidProperties]
    [./simple_fluid]
      type = Water97FluidProperties 
    [../]
  [../]
[]

[Materials]
  [./temperature]
    type = PorousFlowTemperature
    temperature = T
  [../]
  [./temperature_nodal]
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
  [./simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
    at_nodes = true
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
    at_nodes = true
  [../]
  [./simple_fluid_qp]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  [../]
  [./thermal_conductivity_matrix]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.58 0 0 0 2.58 0 0 0 2.58'
    block = 'matrix1 matrix2'
  [../]
  [./thermal_conductivity_fracture]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '0.68 0 0 0 0.68 0 0 0 0.68'
    block = 'fracture'
  [../]
  [./poro_fracture]
    type = PorousFlowPorosityConst
    porosity = 0.05   # = a * phif
    block = 'fracture'
  [../]
  [./poro_matrix]
    type = PorousFlowPorosityConst
    porosity = 0.01
    block = 'matrix1 matrix2'
  [../]
  [./poro_fracture_nodal]
    type = PorousFlowPorosityConst
    at_nodes = true
    porosity = 0.05  # = a * phif
    block = 'fracture'
  [../]
  [./poro_matrix_nodal]
    type = PorousFlowPorosityConst
    at_nodes = true
    porosity = 0.01
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
    type = PorousFlowPermeabilityConst
    permeability = '1.05e-8 0 0 0 1.05e-8 0 0 0 1.05e-8'   # 1.8e-11 = a * kf
    block = 'fracture'
  [../]
  [./permeability_matrix]
    type = PorousFlowPermeabilityConst
    permeability = '1.096e-16 0 0 0 1.096e-16 0 0 0 1.096e-16'
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
  [../]
  [./thermal_expansion_strain]
    type = ComputeThermalExpansionEigenstrain
    stress_free_temperature = 293
    thermal_expansion_coeff = 1.0E-5 
    temperature = T
    eigenstrain_name = thermal_strain 
  [../]

  [./strain_from_initial_stress]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = 'kxx_fcn 0 0
                      0 kyy_fcn 0
                      0 0 weight_fcn'
    eigenstrain_name = ini_stress
  [../]

  [./strain]
    type = ComputeSmallStrain
    displacements = 'sdisp_x sdisp_y sdisp_z'
    eigenstrain_names = 'ini_stress thermal_strain'
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
  [../]
[]

[Preconditioning]
  active = 'superlu'
  [./basic]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2             '
  [../]
  [./superlu]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
    petsc_options_value = 'gmres lu superlu_dist'
  [../]
  [./sub_pc_superlu]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_mat_solver_package'
    petsc_options_value = 'gmres asm lu NONZERO superlu_dist'
  [../]
  [./mumps]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  [../]
[]

[Executioner]
 # type = Transient
 # solve_type = NEWTON
 # num_steps = 1
 # dt = 86400.0
 # 
 # l_max_its = 20
 # l_tol = 1e-5
 # nl_max_its = 30
 # nl_abs_tol = 1e-9
  type = Steady
  solve_type = 'NEWTON'
  l_max_its = 500
  l_tol = 1e-4
  nl_max_its = 300
  nl_abs_tol = 1e-6
[]


[Outputs]
 # [./csv]
 #   type = CSV
 #   execute_on = 'final'
 # [../]
  [./exduos]
    type = Exodus
    file_base = steady_out_fine
    execute_on = 'timestep_end'
  [../]
[]
