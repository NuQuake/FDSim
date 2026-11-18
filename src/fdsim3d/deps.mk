# Module dependencies (generated from the original Makefile; regenerate with e.g. makedepf90 or fortdepend if USE statements change)
$(O)/alloc_nr.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_nr_bound.o $(O)/mod_precision.o
$(O)/alloc_t_seis.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_mpi_shared.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_wavefield.o
$(O)/aux1.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o
$(O)/bcucof.o: $(O)/mod_precision.o
$(O)/coarse2fine_txz.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_precision.o $(O)/mod_stressfield.o
$(O)/coarse2fine_tyz.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_precision.o $(O)/mod_stressfield.o
$(O)/coarse2fine_w.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_interfaces.o \
	$(O)/mod_precision.o $(O)/mod_wavefield.o
$(O)/communicate.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_mpi_shared.o $(O)/mod_precision.o
$(O)/communicate_f.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o
$(O)/communicate_plane.o: $(O)/mod_auxil.o $(O)/mod_mpi_shared.o $(O)/mod_precision.o
$(O)/communicate_stressfield.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/communicate_wavefield.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_precision.o \
	$(O)/mod_wavefield.o
$(O)/compliances.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/distribute.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_mpi_shared.o
$(O)/domain_decomposition.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o
$(O)/energy.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/fdl2mt.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/fine2coarse.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_mpi_shared.o $(O)/mod_precision.o
$(O)/higdce2_fd.o: $(O)/mod_precision.o
$(O)/higdce_fd.o: $(O)/mod_precision.o
$(O)/higdon_fd.o: $(O)/mod_precision.o
$(O)/hplane_0_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_0_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_0f_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_0f_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_1_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_1_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_1f_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_1f_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_ic_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_ic_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_if_stress.o: $(O)/mod_auxil.o $(O)/mod_grid_med.o $(O)/mod_interfaces.o \
	$(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o $(O)/mod_strainfield.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_if_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_lin_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_lin_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_lin_w_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_lin_w_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_linpm1_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_linpm1_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_lpasf_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_interfaces.o \
	$(O)/mod_stressfield.o
$(O)/hplane_lpasm1f_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_lpasm1f_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_lpasm1m2_stress.o: $(O)/mod_contr_data.o $(O)/mod_interfaces.o $(O)/mod_stressfield.o
$(O)/hplane_lpasm1m2_velocity.o: $(O)/mod_contr_data.o $(O)/mod_interfaces.o $(O)/mod_wavefield.o
$(O)/hplane_mz1_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_pml_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_pml_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_pmlzh0_stress.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_interfaces.o $(O)/mod_nr_bound.o $(O)/mod_pml.o $(O)/mod_precision.o \
	$(O)/mod_strainfield.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
$(O)/hplane_pmlzh0_velocity.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_nr_bound.o \
	$(O)/mod_precision.o $(O)/mod_wavefield.o
$(O)/inh_nr.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_nr_bound.o $(O)/mod_precision.o
$(O)/inh_v_u_2nd.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_u_2nd_0.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_u_int_2nd.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_u_int_2nd_pml.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_uv_int.o: $(O)/mod_auxil.o $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_uv_int_0.o: $(O)/mod_auxil.o $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_uv_int_1_si.o: $(O)/mod_precision.o
$(O)/inh_v_v_2nd.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_v_2nd_0.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_v_int_2nd.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_v_int_2nd_pml.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_w_2nd.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_w_2nd_0.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_w_int.o: $(O)/mod_auxil.o $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_w_int_0.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_w_int_1.o: $(O)/mod_auxil.o $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_w_int_1_si.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_w_int_2nd.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/inh_v_w_int_2nd_pml.o: $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/interpol.o: $(O)/mod_precision.o
$(O)/lanczos.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_precision.o
$(O)/main.o: $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_interfaces.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o $(O)/mod_time.o
$(O)/mod_auxil.o: $(O)/mod_nr_bound.o $(O)/mod_precision.o
$(O)/mod_contr_data.o: $(O)/mod_precision.o
$(O)/mod_grid_med.o: $(O)/mod_precision.o
$(O)/mod_interfaces.o: $(O)/mod_auxil.o $(O)/mod_grid_med.o $(O)/mod_nr_bound.o $(O)/mod_precision.o
$(O)/mod_nr_bound.o: $(O)/mod_precision.o
$(O)/mod_pml.o: $(O)/mod_precision.o
$(O)/mod_strainfield.o: $(O)/mod_precision.o
$(O)/mod_stressfield.o: $(O)/mod_precision.o
$(O)/mod_wavefield.o: $(O)/mod_precision.o
$(O)/nonref.o: $(O)/mod_contr_data.o $(O)/mod_precision.o
$(O)/open_sr_file.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o
$(O)/OptAbsBc.o: $(O)/mod_precision.o
$(O)/ppml.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_pml.o $(O)/mod_precision.o
$(O)/ppmlf.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_pml.o \
	$(O)/mod_precision.o
$(O)/pre_nr_coef1.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_nr_bound.o \
	$(O)/mod_precision.o
$(O)/pre_nr_coef2.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_nr_bound.o \
	$(O)/mod_precision.o
$(O)/read_model.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o
$(O)/read_model_par.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o \
	$(O)/mod_mpi_shared.o $(O)/mod_precision.o
$(O)/read_open.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o
$(O)/src_add.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_precision.o \
	$(O)/mod_wavefield.o
$(O)/src_v.o: $(O)/mod_auxil.o $(O)/mod_grid_med.o $(O)/mod_precision.o
$(O)/sto_bor.o: $(O)/mod_nr_bound.o $(O)/mod_precision.o
$(O)/store_seismograms.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o $(O)/mod_wavefield.o
$(O)/store_snapshots.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o $(O)/mod_wavefield.o
$(O)/store_snapshots_f.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_mpi_shared.o \
	$(O)/mod_precision.o $(O)/mod_wavefield.o
$(O)/strain_xx.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_xy.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_xz.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_xz_1.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_xz_2nd.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_xz_2nd_pml.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_yy.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_yz.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_yz_1.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_yz_2nd.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_yz_2nd_pml.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_zz.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_zz_0.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_zz_2nd.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/strain_zz_2nd_pml.o: $(O)/mod_auxil.o $(O)/mod_interfaces.o $(O)/mod_precision.o
$(O)/stress_nn.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_interfaces.o \
	$(O)/mod_precision.o
$(O)/stress_xy.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_interfaces.o \
	$(O)/mod_precision.o
$(O)/stress_xz.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_interfaces.o \
	$(O)/mod_precision.o
$(O)/stress_yz.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_interfaces.o \
	$(O)/mod_precision.o
$(O)/time_loop.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_interfaces.o \
	$(O)/mod_mpi_shared.o $(O)/mod_nr_bound.o $(O)/mod_t_loop.o $(O)/mod_time.o \
	$(O)/mod_wavefield.o
$(O)/tl_alloc.o: $(O)/mod_auxil.o $(O)/mod_contr_data.o $(O)/mod_grid_med.o $(O)/mod_stressfield.o \
	$(O)/mod_wavefield.o
$(O)/tl_allocf.o: $(O)/mod_auxil.o $(O)/mod_stressfield.o $(O)/mod_wavefield.o
