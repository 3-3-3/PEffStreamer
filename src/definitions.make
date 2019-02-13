OBJS := m_units_constants.o m_config.o m_lookup_table.o m_random.o		\
	m_photoi_mc.o m_streamer.o m_geometry.o m_transport_data.o m_field.o	\
	m_init_cond.o m_photoi_helmh.o m_photoi.o m_chemistry.o m_types.o	\
	m_gas.o m_refine.o m_fluid_lfa.o m_dt.o m_advance.o m_advance_base.o	\
	m_user_methods.o m_table_data.o m_output.o

# Hide some incorrect warnings
m_photoi_helmh.o: FFLAGS += -Wno-unused-function
m_photoi.o: FFLAGS += -Wno-unused-function

# Dependency information
m_advance_base.o: m_config.mod
m_advance_base.o: m_types.mod
m_advance.o: m_advance_base.mod
m_advance.o: m_chemistry.mod
m_advance.o: m_dt.mod
m_advance.o: m_field.mod
m_advance.o: m_fluid_lfa.mod
m_advance.o: m_streamer.mod
m_chemistry.o: m_advance_base.mod
m_chemistry.o: m_config.mod
m_chemistry.o: m_gas.mod
m_chemistry.o: m_lookup_table.mod
m_chemistry.o: m_table_data.mod
m_chemistry.o: m_transport_data.mod
m_chemistry.o: m_types.mod
m_chemistry.o: m_units_constants.mod
m_dt.o: m_config.mod
m_field.o: m_chemistry.mod
m_field.o: m_config.mod
m_field.o: m_lookup_table.mod
m_field.o: m_streamer.mod
m_field.o: m_table_data.mod
m_field.o: m_types.mod
m_field.o: m_units_constants.mod
m_fluid_lfa.o: m_chemistry.mod
m_fluid_lfa.o: m_dt.mod
m_fluid_lfa.o: m_gas.mod
m_fluid_lfa.o: m_lookup_table.mod
m_fluid_lfa.o: m_photoi.mod
m_fluid_lfa.o: m_streamer.mod
m_fluid_lfa.o: m_transport_data.mod
m_fluid_lfa.o: m_units_constants.mod
m_gas.o: m_config.mod
m_gas.o: m_types.mod
m_gas.o: m_units_constants.mod
m_init_cond.o: m_chemistry.mod
m_init_cond.o: m_config.mod
m_init_cond.o: m_gas.mod
m_init_cond.o: m_geometry.mod
m_init_cond.o: m_streamer.mod
m_init_cond.o: m_types.mod
m_output.o: m_advance.mod
m_output.o: m_config.mod
m_output.o: m_field.mod
m_output.o: m_photoi.mod
m_output.o: m_streamer.mod
m_output.o: m_types.mod
m_output.o: m_user_methods.mod
m_photoi_helmh.o: m_config.mod
m_photoi_helmh.o: m_gas.mod
m_photoi_helmh.o: m_streamer.mod
m_photoi_helmh.o: m_units_constants.mod
m_photoi_mc.o: m_config.mod
m_photoi_mc.o: m_gas.mod
m_photoi_mc.o: m_lookup_table.mod
m_photoi_mc.o: m_random.mod
m_photoi_mc.o: m_streamer.mod
m_photoi_mc.o: m_units_constants.mod
m_photoi.o: m_config.mod
m_photoi.o: m_gas.mod
m_photoi.o: m_lookup_table.mod
m_photoi.o: m_photoi_helmh.mod
m_photoi.o: m_photoi_mc.mod
m_photoi.o: m_streamer.mod
m_photoi.o: m_transport_data.mod
m_photoi.o: m_types.mod
m_photoi.o: m_units_constants.mod
m_refine.o: m_config.mod
m_refine.o: m_gas.mod
m_refine.o: m_geometry.mod
m_refine.o: m_init_cond.mod
m_refine.o: m_lookup_table.mod
m_refine.o: m_streamer.mod
m_refine.o: m_transport_data.mod
m_streamer.o: m_chemistry.mod
m_streamer.o: m_config.mod
m_streamer.o: m_gas.mod
m_streamer.o: m_lookup_table.mod
m_streamer.o: m_random.mod
m_streamer.o: m_types.mod
m_streamer.o: m_units_constants.mod
m_table_data.o: m_config.mod
m_table_data.o: m_types.mod
m_transport_data.o: m_config.mod
m_transport_data.o: m_gas.mod
m_transport_data.o: m_lookup_table.mod
m_transport_data.o: m_table_data.mod
m_transport_data.o: m_types.mod
streamer.o: m_advance_base.mod
streamer.o: m_advance.mod
streamer.o: m_chemistry.mod
streamer.o: m_config.mod
streamer.o: m_dt.mod
streamer.o: m_field.mod
streamer.o: m_fluid_lfa.mod
streamer.o: m_gas.mod
streamer.o: m_init_cond.mod
streamer.o: m_output.mod
streamer.o: m_photoi.mod
streamer.o: m_refine.mod
streamer.o: m_streamer.mod
streamer.o: m_table_data.mod
streamer.o: m_transport_data.mod
streamer.o: m_types.mod
streamer.o: m_user_methods.mod
streamer.o: m_user.mod
