# Module dependencies (generated from the original Makefile; regenerate with e.g. makedepf90 or fortdepend if USE statements change)
$(O)/add_new_cell_type.o: $(O)/mod_model.o $(O)/mod_parameterization.o
$(O)/allocate_arrays.o: $(O)/mod_const.o $(O)/mod_model.o $(O)/mod_parameterization.o
$(O)/enlarge.o: $(O)/mod_const.o $(O)/mod_model.o $(O)/mod_parameterization.o $(O)/nrtype.o
$(O)/fdmod3d.o: $(O)/mod_interfaces.o $(O)/mod_model.o $(O)/mod_mpi_shared.o
$(O)/filter.o: $(O)/mod_model.o $(O)/nrtype.o
$(O)/find_max_depth.o: $(O)/mod_model.o $(O)/nrtype.o
$(O)/func_mu_la_rho.o: $(O)/mod_func.o $(O)/mod_interfaces.o $(O)/mod_model.o \
	$(O)/mod_parameterization.o $(O)/nrtype.o
$(O)/func_ymu_yl2m.o: $(O)/mod_atten.o $(O)/mod_const.o $(O)/mod_func.o $(O)/mod_interfaces.o \
	$(O)/mod_model.o $(O)/mod_parameterization.o $(O)/nrtype.o
$(O)/get_fd_cell_index.o: $(O)/mod_interfaces.o $(O)/mod_model.o $(O)/mod_parameterization.o
$(O)/get_layer.o: $(O)/mod_interfaces.o $(O)/mod_model.o $(O)/nrtype.o
$(O)/getdepth.o: $(O)/mod_model.o $(O)/nrtype.o
$(O)/int3d.o: $(O)/mod_integration.o $(O)/nrtype.o
$(O)/is_in_1D_medium.o: $(O)/mod_interfaces.o $(O)/mod_model.o $(O)/mod_parameterization.o \
	$(O)/nrtype.o
$(O)/mod_atten.o: $(O)/mod_const.o $(O)/nrtype.o
$(O)/mod_const.o: $(O)/nrtype.o
$(O)/mod_effmatpar.o: $(O)/mod_integration.o $(O)/nrtype.o
$(O)/mod_effmatpar_y.o: $(O)/mod_const.o $(O)/mod_integration.o $(O)/nrtype.o
$(O)/mod_func.o: $(O)/mod_interfaces.o $(O)/mod_model.o $(O)/nrtype.o
$(O)/mod_integration.o: $(O)/nrtype.o
$(O)/mod_interfaces.o: $(O)/mod_const.o $(O)/nrtype.o
$(O)/mod_model.o: $(O)/nrtype.o
$(O)/mod_nr.o: $(O)/nrtype.o
$(O)/mod_parameterization.o: $(O)/mod_const.o $(O)/nrtype.o
$(O)/model_lin.o: $(O)/mod_files.o $(O)/mod_model.o $(O)/mod_mpi_shared.o $(O)/nrtype.o
$(O)/par_ha.o: $(O)/mod_const.o $(O)/mod_files.o $(O)/mod_model.o $(O)/mod_parameterization.o
$(O)/parameterization.o: $(O)/mod_interfaces.o $(O)/mod_model.o $(O)/mod_mpi_shared.o \
	$(O)/mod_parameterization.o
$(O)/parameterize_coarser_grid.o: $(O)/mod_atten.o $(O)/mod_const.o $(O)/mod_effmatpar.o \
	$(O)/mod_effmatpar_y.o $(O)/mod_integration.o $(O)/mod_interfaces.o $(O)/mod_model.o \
	$(O)/mod_mpi_shared.o $(O)/mod_parameterization.o $(O)/nrtype.o
$(O)/parameterize_finer_grid.o: $(O)/mod_atten.o $(O)/mod_const.o $(O)/mod_effmatpar.o \
	$(O)/mod_effmatpar_y.o $(O)/mod_integration.o $(O)/mod_interfaces.o $(O)/mod_model.o \
	$(O)/mod_mpi_shared.o $(O)/mod_parameterization.o $(O)/nrtype.o
$(O)/prepare_fine_par.o: $(O)/mod_integration.o $(O)/mod_interfaces.o $(O)/mod_model.o $(O)/nrtype.o
$(O)/prepare_fine_par_y.o: $(O)/mod_atten.o $(O)/mod_const.o $(O)/mod_integration.o \
	$(O)/mod_interfaces.o $(O)/mod_model.o $(O)/mod_parameterization.o $(O)/nrtype.o
$(O)/prepare_input.o: $(O)/mod_files.o $(O)/mod_model.o $(O)/mod_parameterization.o $(O)/nrtype.o
$(O)/prepare_int.o: $(O)/mod_const.o $(O)/mod_integration.o $(O)/mod_model.o $(O)/nrtype.o
$(O)/prepare_y.o: $(O)/mod_atten.o $(O)/mod_const.o $(O)/mod_files.o $(O)/mod_parameterization.o \
	$(O)/nrtype.o
$(O)/print_info.o: $(O)/mod_model.o $(O)/mod_parameterization.o
$(O)/qylsm.o: $(O)/mod_atten.o $(O)/mod_const.o $(O)/nrtype.o
$(O)/read_GBM3D.o: $(O)/mod_model.o $(O)/mod_mpi_shared.o
$(O)/setparams.o: $(O)/mod_const.o $(O)/mod_files.o $(O)/mod_integration.o $(O)/mod_interfaces.o \
	$(O)/mod_model.o $(O)/mod_mpi_shared.o $(O)/mod_parameterization.o $(O)/nrtype.o
$(O)/test_cells.o: $(O)/mod_const.o $(O)/mod_parameterization.o $(O)/nrtype.o
$(O)/test_cellsy.o: $(O)/mod_const.o $(O)/mod_parameterization.o $(O)/nrtype.o
$(O)/test_parameters.o: $(O)/mod_files.o $(O)/mod_integration.o $(O)/mod_model.o $(O)/nrtype.o
$(O)/write_grd.o: $(O)/mod_files.o $(O)/mod_model.o $(O)/nrtype.o
$(O)/write_param.o: $(O)/mod_const.o $(O)/mod_files.o $(O)/mod_model.o $(O)/mod_mpi_shared.o \
	$(O)/mod_parameterization.o $(O)/nrtype.o
