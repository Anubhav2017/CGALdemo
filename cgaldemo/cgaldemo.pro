QT += core
QT -= gui

CONFIG += c++11

TARGET = cgaldemo
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
#    CGAL/CORE/poly/Curves.tcc \
#    CGAL/CORE/poly/Poly.tcc

FORMS += \
    CGAL/Qt/resources/ImageInterface.ui

DISTFILES += \
    CGAL/Qt/resources/qglviewer-icon.xpm \
    CGAL/license/gpl_package_list.txt

HEADERS += \
#    CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h \
#    CGAL/Algebraic_kernel_d/algebraic_curve_kernel_2_tools.h \
#    CGAL/Algebraic_kernel_d/Algebraic_real_d_1.h \
#    CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h \
#    CGAL/Algebraic_kernel_d/Algebraic_real_rep.h \
#    CGAL/Algebraic_kernel_d/Algebraic_real_rep_bfi.h \
#    CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h \
#    CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel_at_alpha.h \
#    CGAL/Algebraic_kernel_d/Bitstream_descartes.h \
#    CGAL/Algebraic_kernel_d/Bitstream_descartes_E08_tree.h \
#    CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree.h \
#    CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h \
#    CGAL/Algebraic_kernel_d/bound_between_1.h \
#    CGAL/Algebraic_kernel_d/construct_binary.h \
#    CGAL/Algebraic_kernel_d/Curve_analysis_2.h \
#    CGAL/Algebraic_kernel_d/Curve_pair_analysis_2.h \
#    CGAL/Algebraic_kernel_d/Descartes.h \
#    CGAL/Algebraic_kernel_d/enums.h \
#    CGAL/Algebraic_kernel_d/Event_line_builder.h \
#    CGAL/Algebraic_kernel_d/exceptions.h \
#    CGAL/Algebraic_kernel_d/flags.h \
#    CGAL/Algebraic_kernel_d/Float_traits.h \
#    CGAL/Algebraic_kernel_d/Interval_evaluate_1.h \
#    CGAL/Algebraic_kernel_d/Interval_evaluate_2.h \
#    CGAL/Algebraic_kernel_d/LRU_hashed_map.h \
#    CGAL/Algebraic_kernel_d/macros.h \
#    CGAL/Algebraic_kernel_d/Real_embeddable_extension.h \
#    CGAL/Algebraic_kernel_d/Real_roots.h \
#    CGAL/Algebraic_kernel_d/refine_zero_against.h \
#    CGAL/Algebraic_kernel_d/shear.h \
#    CGAL/Algebraic_kernel_d/Shear_controller.h \
#    CGAL/Algebraic_kernel_d/Shear_transformation.h \
#    CGAL/Algebraic_kernel_d/Status_line_CA_1.h \
#    CGAL/Algebraic_kernel_d/Status_line_CPA_1.h \
#    CGAL/Algebraic_kernel_d/univariate_polynomial_utils.h \
#    CGAL/Algebraic_kernel_d/Xy_coordinate_2.h \
#    CGAL/Algebraic_kernel_for_circles/function_objects_on_roots_and_polynomials_2_2.h \
#    CGAL/Algebraic_kernel_for_circles/internal_functions_comparison_root_for_circles_2_2.h \
#    CGAL/Algebraic_kernel_for_circles/internal_functions_on_roots_and_polynomial_1_2_and_2_2.h \
#    CGAL/Algebraic_kernel_for_circles/internal_functions_on_roots_and_polynomials_2_2.h \
#    CGAL/Algebraic_kernel_for_spheres/function_objects_on_roots_and_polynomials_2_3.h \
#    CGAL/Algebraic_kernel_for_spheres/internal_functions_comparison_root_for_spheres_2_3.h \
#    CGAL/Algebraic_kernel_for_spheres/internal_functions_on_roots_and_polynomial_1_3_and_2_3.h \
#    CGAL/Algebraic_kernel_for_spheres/internal_functions_on_roots_and_polynomials_1_3.h \
#    CGAL/Algebraic_kernel_for_spheres/internal_functions_on_roots_and_polynomials_2_3.h \
#    CGAL/Apollonius_graph_2/uncertain/uncertain_functions_on_signs.h \
#    CGAL/Apollonius_graph_2/uncertain/Uncertain_is_hidden_C2.h \
#    CGAL/Apollonius_graph_2/uncertain/Uncertain_oriented_side_of_bisector_C2.h \
#    CGAL/Apollonius_graph_2/uncertain/Uncertain_vertex_conflict_2.h \
#    CGAL/Apollonius_graph_2/Apollonius_graph_2_impl.h \
#    CGAL/Apollonius_graph_2/Apollonius_graph_hierarchy_2_impl.h \
#    CGAL/Apollonius_graph_2/basic.h \
#    CGAL/Apollonius_graph_2/Bounded_side_of_ccw_circle_C2.h \
#    CGAL/Apollonius_graph_2/check_filter.h \
#    CGAL/Apollonius_graph_2/comparator_profiler.h \
#    CGAL/Apollonius_graph_2/compare_quadratic.h \
#    CGAL/Apollonius_graph_2/Compare_weight_2.h \
#    CGAL/Apollonius_graph_2/Compare_x_2.h \
#    CGAL/Apollonius_graph_2/Compare_y_2.h \
#    CGAL/Apollonius_graph_2/Constructions_C2.h \
#    CGAL/Apollonius_graph_2/Constructions_ftC2.h \
#    CGAL/Apollonius_graph_2/Constructions_rtH2.h \
#    CGAL/Apollonius_graph_2/Finite_edge_test8_C2.h \
#    CGAL/Apollonius_graph_2/Finite_edge_test_C2.h \
#    CGAL/Apollonius_graph_2/Incircle8_C2.h \
#    CGAL/Apollonius_graph_2/Incircle_C2.h \
#    CGAL/Apollonius_graph_2/Infinite_edge_test_C2.h \
#    CGAL/Apollonius_graph_2/Is_degenerate_edge_C2.h \
#    CGAL/Apollonius_graph_2/Is_hidden_C2.h \
#    CGAL/Apollonius_graph_2/Kernel_wrapper_2.h \
#    CGAL/Apollonius_graph_2/Orientation8_C2.h \
#    CGAL/Apollonius_graph_2/Orientation_2.h \
#    CGAL/Apollonius_graph_2/Oriented_side_of_bisector_C2.h \
#    CGAL/Apollonius_graph_2/Predicate_constructions_C2.h \
#    CGAL/Apollonius_graph_2/predicate_profiler.h \
#    CGAL/Apollonius_graph_2/Predicates_C2.h \
#    CGAL/Apollonius_graph_2/Traits_wrapper_2.h \
#    CGAL/Approximate_min_ellipsoid_d/Approximate_min_ellipsoid_d_configure.h \
#    CGAL/Approximate_min_ellipsoid_d/Approximate_min_ellipsoid_d_debug.h \
#    CGAL/Approximate_min_ellipsoid_d/Approximate_min_ellipsoid_d_impl.h \
#    CGAL/Approximate_min_ellipsoid_d/Khachiyan_approximation.h \
#    CGAL/Approximate_min_ellipsoid_d/Khachiyan_approximation_impl.h \
#    CGAL/Arithmetic_kernel/Arithmetic_kernel_base.h \
#    CGAL/Arr_geometry_traits/IO/Polycurve_2_iostream.h \
#    CGAL/Arr_geometry_traits/Arr_plane_3.h \
#    CGAL/Arr_geometry_traits/Bezier_bounding_rational_traits.h \
#    CGAL/Arr_geometry_traits/Bezier_cache.h \
#    CGAL/Arr_geometry_traits/Bezier_curve_2.h \
#    CGAL/Arr_geometry_traits/Bezier_point_2.h \
#    CGAL/Arr_geometry_traits/Bezier_x_monotone_2.h \
#    CGAL/Arr_geometry_traits/Circle_segment_2.h \
#    CGAL/Arr_geometry_traits/Conic_arc_2.h \
#    CGAL/Arr_geometry_traits/Conic_intersections_2.h \
#    CGAL/Arr_geometry_traits/Conic_point_2.h \
#    CGAL/Arr_geometry_traits/Conic_x_monotone_arc_2.h \
#    CGAL/Arr_geometry_traits/Consolidated_curve_data_aux.h \
#    CGAL/Arr_geometry_traits/Curve_data_aux.h \
#    CGAL/Arr_geometry_traits/de_Casteljau_2.h \
#    CGAL/Arr_geometry_traits/One_root_number.h \
#    CGAL/Arr_geometry_traits/Polycurve_2.h \
#    CGAL/Arr_geometry_traits/Polyline_2.h \
#    CGAL/Arr_geometry_traits/Rational_arc_2.h \
#    CGAL/Arr_geometry_traits/Segment_assertions.h \
#    CGAL/Arr_point_location/Arr_batched_point_location_traits_2.h \
#    CGAL/Arr_point_location/Arr_landmarks_pl_impl.h \
#    CGAL/Arr_point_location/Arr_lm_generator_base.h \
#    CGAL/Arr_point_location/Arr_lm_grid_generator.h \
#    CGAL/Arr_point_location/Arr_lm_halton_generator.h \
#    CGAL/Arr_point_location/Arr_lm_middle_edges_generator.h \
#    CGAL/Arr_point_location/Arr_lm_nearest_neighbor.h \
#    CGAL/Arr_point_location/Arr_lm_random_generator.h \
#    CGAL/Arr_point_location/Arr_lm_specified_points_generator.h \
#    CGAL/Arr_point_location/Arr_lm_vertices_generator.h \
#    CGAL/Arr_point_location/Arr_naive_point_location_impl.h \
#    CGAL/Arr_point_location/Arr_simple_point_location_impl.h \
#    CGAL/Arr_point_location/Arr_trapezoid_ric_pl_impl.h \
#    CGAL/Arr_point_location/Arr_triangulation_pl_functions.h \
#    CGAL/Arr_point_location/Arr_triangulation_pl_impl.h \
#    CGAL/Arr_point_location/Arr_walk_along_line_pl_impl.h \
#    CGAL/Arr_point_location/Td_active_edge.h \
#    CGAL/Arr_point_location/Td_active_fictitious_vertex.h \
#    CGAL/Arr_point_location/Td_active_trapezoid.h \
#    CGAL/Arr_point_location/Td_active_vertex.h \
#    CGAL/Arr_point_location/Td_dag.h \
#    CGAL/Arr_point_location/Td_dag_node.h \
#    CGAL/Arr_point_location/Td_inactive_edge.h \
#    CGAL/Arr_point_location/Td_inactive_fictitious_vertex.h \
#    CGAL/Arr_point_location/Td_inactive_trapezoid.h \
#    CGAL/Arr_point_location/Td_inactive_vertex.h \
#    CGAL/Arr_point_location/Td_ninetuple.h \
#    CGAL/Arr_point_location/Td_predicates.h \
#    CGAL/Arr_point_location/Td_traits.h \
#    CGAL/Arr_point_location/Td_X_trapezoid.h \
#    CGAL/Arr_point_location/Trapezoidal_decomposition_2.h \
#    CGAL/Arr_point_location/Trapezoidal_decomposition_2_impl.h \
#    CGAL/Arr_point_location/Trapezoidal_decomposition_2_iostream.h \
#    CGAL/Arr_point_location/Trapezoidal_decomposition_2_misc.h \
#    CGAL/Arr_rat_arc/Algebraic_point_2.h \
#    CGAL/Arr_rat_arc/Base_rational_arc_ds_1.h \
#    CGAL/Arr_rat_arc/Cache.h \
#    CGAL/Arr_rat_arc/Rational_arc_d_1.h \
#    CGAL/Arr_rat_arc/Rational_function.h \
#    CGAL/Arr_rat_arc/Rational_function_canonicalized_pair.h \
#    CGAL/Arr_rat_arc/Rational_function_ordered_pair.h \
#    CGAL/Arr_rat_arc/Rational_function_pair.h \
#    CGAL/Arr_rat_arc/Singleton.h \
#    CGAL/Arr_spherical_gaussian_map_3/Arr_on_sphere_transformation.h \
#    CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm.h \
#    CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_arr_dcel.h \
#    CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_initializer_visitor.h \
#    CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_overlay.h \
#    CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_polyhedron_3.h \
#    CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_traits.h \
#    CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_transformation.h \
#    CGAL/Arr_spherical_gaussian_map_3/Arr_spherical_gaussian_map_3.h \
#    CGAL/Arr_spherical_gaussian_map_3/Arr_transform_on_sphere.h \
#    CGAL/Arr_topology_traits/Arr_bounded_planar_batched_pl_helper.h \
#    CGAL/Arr_topology_traits/Arr_bounded_planar_construction_helper.h \
#    CGAL/Arr_topology_traits/Arr_bounded_planar_insertion_helper.h \
#    CGAL/Arr_topology_traits/Arr_bounded_planar_overlay_helper.h \
#    CGAL/Arr_topology_traits/Arr_bounded_planar_topology_traits_2_impl.h \
#    CGAL/Arr_topology_traits/Arr_bounded_planar_vert_decomp_helper.h \
#    CGAL/Arr_topology_traits/Arr_inc_insertion_zone_visitor.h \
#    CGAL/Arr_topology_traits/Arr_planar_topology_traits_base_2.h \
#    CGAL/Arr_topology_traits/Arr_spherical_batched_pl_helper.h \
#    CGAL/Arr_topology_traits/Arr_spherical_construction_helper.h \
#    CGAL/Arr_topology_traits/Arr_spherical_insertion_helper.h \
#    CGAL/Arr_topology_traits/Arr_spherical_overlay_helper.h \
#    CGAL/Arr_topology_traits/Arr_spherical_topology_traits_2_impl.h \
#    CGAL/Arr_topology_traits/Arr_spherical_vert_decomp_helper.h \
#    CGAL/Arr_topology_traits/Arr_unb_planar_batched_pl_helper.h \
#    CGAL/Arr_topology_traits/Arr_unb_planar_construction_helper.h \
#    CGAL/Arr_topology_traits/Arr_unb_planar_insertion_helper.h \
#    CGAL/Arr_topology_traits/Arr_unb_planar_overlay_helper.h \
#    CGAL/Arr_topology_traits/Arr_unb_planar_topology_traits_2_impl.h \
#    CGAL/Arr_topology_traits/Arr_unb_planar_vert_decomp_helper.h \
#    CGAL/Arrangement_2/Arr_compute_zone_visitor.h \
#    CGAL/Arrangement_2/Arr_default_planar_topology.h \
#    CGAL/Arrangement_2/Arr_do_intersect_zone_visitor.h \
#    CGAL/Arrangement_2/Arr_on_surface_with_history_2_impl.h \
#    CGAL/Arrangement_2/Arr_traits_adaptor_2.h \
#    CGAL/Arrangement_2/Arr_traits_adaptor_2_dispatching.h \
#    CGAL/Arrangement_2/Arr_with_history_accessor.h \
#    CGAL/Arrangement_2/Arrangement_2_iterators.h \
#    CGAL/Arrangement_2/Arrangement_on_surface_2_global.h \
#    CGAL/Arrangement_2/Arrangement_on_surface_2_impl.h \
#    CGAL/Arrangement_2/arrangement_type_traits.h \
#    CGAL/Arrangement_2/Arrangement_zone_2_impl.h \
#    CGAL/Arrangement_2/graph_traits_dual.h \
#    CGAL/auto_link/auto_link.h \
#    CGAL/auto_link/CGAL.h \
#    CGAL/auto_link/CORE.h \
#    CGAL/auto_link/ImageIO.h \
#    CGAL/auto_link/LAPACK.h \
#    CGAL/auto_link/Qt.h \
#    CGAL/auto_link/TAUCS.h \
#    CGAL/Barycentric_coordinates_2/barycentric_enum_2.h \
#    CGAL/Barycentric_coordinates_2/Discrete_harmonic_2.h \
#    CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h \
#    CGAL/Barycentric_coordinates_2/Mean_value_2.h \
#    CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h \
#    CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h \
#    CGAL/Barycentric_coordinates_2/Wachspress_2.h \
#    CGAL/Boolean_set_operations_2/Bso_internal_functions.h \
#    CGAL/Boolean_set_operations_2/Ccb_curve_iterator.h \
#    CGAL/Boolean_set_operations_2/Curve_with_halfedge.h \
#    CGAL/Boolean_set_operations_2/Gps_agg_meta_traits.h \
#    CGAL/Boolean_set_operations_2/Gps_agg_op.h \
#    CGAL/Boolean_set_operations_2/Gps_agg_op_surface_sweep_2.h \
#    CGAL/Boolean_set_operations_2/Gps_agg_op_visitor.h \
#    CGAL/Boolean_set_operations_2/Gps_base_functor.h \
#    CGAL/Boolean_set_operations_2/Gps_bfs_base_visitor.h \
#    CGAL/Boolean_set_operations_2/Gps_bfs_intersection_visitor.h \
#    CGAL/Boolean_set_operations_2/Gps_bfs_join_visitor.h \
#    CGAL/Boolean_set_operations_2/Gps_bfs_scanner.h \
#    CGAL/Boolean_set_operations_2/Gps_bfs_xor_visitor.h \
#    CGAL/Boolean_set_operations_2/Gps_default_dcel.h \
#    CGAL/Boolean_set_operations_2/Gps_default_traits.h \
#    CGAL/Boolean_set_operations_2/Gps_difference_functor.h \
#    CGAL/Boolean_set_operations_2/Gps_do_intersect_functor.h \
#    CGAL/Boolean_set_operations_2/Gps_insertion_meta_traits.h \
#    CGAL/Boolean_set_operations_2/Gps_intersection_functor.h \
#    CGAL/Boolean_set_operations_2/Gps_join_functor.h \
#    CGAL/Boolean_set_operations_2/Gps_merge.h \
#    CGAL/Boolean_set_operations_2/Gps_on_surface_base_2.h \
#    CGAL/Boolean_set_operations_2/Gps_on_surface_base_2_impl.h \
#    CGAL/Boolean_set_operations_2/Gps_polygon_simplifier.h \
#    CGAL/Boolean_set_operations_2/Gps_polygon_validation.h \
#    CGAL/Boolean_set_operations_2/Gps_simplifier_traits.h \
#    CGAL/Boolean_set_operations_2/Gps_sym_diff_functor.h \
#    CGAL/Boolean_set_operations_2/Gps_traits_adaptor.h \
#    CGAL/Boolean_set_operations_2/Gps_traits_decorator.h \
#    CGAL/Boolean_set_operations_2/Indexed_event.h \
#    CGAL/Boolean_set_operations_2/Point_with_vertex.h \
#    CGAL/Boolean_set_operations_2/Polygon_2_curve_iterator.h \
#    CGAL/boost/bimap/multiset_of.hpp \
#    CGAL/boost/graph/internal/Has_member_clear.h \
#    CGAL/boost/graph/internal/Has_member_id.h \
#    CGAL/boost/graph/internal/helpers.h \
#    CGAL/boost/graph/internal/OM_iterator_from_circulator.h \
#    CGAL/boost/graph/METIS/partition_dual_graph.h \
#    CGAL/boost/graph/METIS/partition_graph.h \
#    CGAL/boost/graph/backward_compatibility_functions.h \
#    CGAL/boost/graph/boost_parameters_interface.h \
#    CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h \
#    CGAL/boost/graph/copy_face_graph.h \
#    CGAL/boost/graph/dijkstra_shortest_paths.h \
#    CGAL/boost/graph/dijkstra_shortest_paths.hpp \
#    CGAL/boost/graph/Dual.h \
#    CGAL/boost/graph/Euler_operations.h \
#    CGAL/boost/graph/Face_filtered_graph.h \
#    CGAL/boost/graph/graph_concepts.h \
#    CGAL/boost/graph/graph_traits_Arrangement_2.h \
#    CGAL/boost/graph/graph_traits_Constrained_Delaunay_triangulation_2.h \
#    CGAL/boost/graph/graph_traits_Constrained_triangulation_2.h \
#    CGAL/boost/graph/graph_traits_Constrained_triangulation_plus_2.h \
#    CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h \
#    CGAL/boost/graph/graph_traits_Dual_Arrangement_2.h \
#    CGAL/boost/graph/graph_traits_HalfedgeDS.h \
#    CGAL/boost/graph/graph_traits_HalfedgeDS_default.h \
#    CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h \
#    CGAL/boost/graph/graph_traits_OpenMesh.h \
#    CGAL/boost/graph/graph_traits_Polyhedron_3.h \
#    CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h \
#    CGAL/boost/graph/graph_traits_Regular_triangulation_2.h \
#    CGAL/boost/graph/graph_traits_Seam_mesh.h \
#    CGAL/boost/graph/graph_traits_Surface_mesh.h \
#    CGAL/boost/graph/graph_traits_Triangulation_2.h \
#    CGAL/boost/graph/graph_traits_Triangulation_data_structure_2.h \
#    CGAL/boost/graph/graph_traits_Triangulation_hierarchy_2.h \
#    CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h \
#    CGAL/boost/graph/Graph_with_descriptor_with_graph.h \
#    CGAL/boost/graph/Graph_with_descriptor_with_graph_fwd.h \
#    CGAL/boost/graph/halfedge_graph_traits.h \
#    CGAL/boost/graph/halfedge_graph_traits_HalfedgeDS.h \
#    CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h \
#    CGAL/boost/graph/helpers.h \
#    CGAL/boost/graph/io.h \
#    CGAL/boost/graph/iterator.h \
#    CGAL/boost/graph/named_function_params.h \
#    CGAL/boost/graph/named_params_helper.h \
#    CGAL/boost/graph/parameters_interface.h \
#    CGAL/boost/graph/partition.h \
#    CGAL/boost/graph/properties.h \
#    CGAL/boost/graph/properties_Linear_cell_complex_for_combinatorial_map.h \
#    CGAL/boost/graph/properties_OpenMesh.h \
#    CGAL/boost/graph/properties_Polyhedron_3.h \
#    CGAL/boost/graph/properties_Polyhedron_3_features.h \
#    CGAL/boost/graph/properties_Polyhedron_3_time_stamp.h \
#    CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h \
#    CGAL/boost/graph/properties_Seam_mesh.h \
#    CGAL/boost/graph/properties_Surface_mesh.h \
#    CGAL/boost/graph/properties_Surface_mesh_features.h \
#    CGAL/boost/graph/properties_Surface_mesh_time_stamp.h \
#    CGAL/boost/graph/properties_TriMesh_ArrayKernelT.h \
#    CGAL/boost/graph/property_maps.h \
#    CGAL/boost/graph/Seam_mesh.h \
#    CGAL/boost/graph/selection.h \
#    CGAL/boost/graph/split_graph_into_polylines.h \
#    CGAL/boost/graph/visitor.h \
#    CGAL/boost/iterator/counting_iterator.hpp \
#    CGAL/boost/iterator/transform_iterator.hpp \
#    CGAL/boost/bimap.hpp \
#    CGAL/boost/parameter.h \
#    CGAL/Box_intersection_d/Box_d.h \
#    CGAL/Box_intersection_d/box_limits.h \
#    CGAL/Box_intersection_d/Box_traits_d.h \
#    CGAL/Box_intersection_d/Box_with_handle_d.h \
#    CGAL/Box_intersection_d/Box_with_info_d.h \
#    CGAL/Box_intersection_d/segment_tree.h \
    CGAL/Cartesian/Aff_transformation_2.h \
    CGAL/Cartesian/Aff_transformation_3.h \
    CGAL/Cartesian/Aff_transformation_rep_2.h \
    CGAL/Cartesian/Aff_transformation_rep_3.h \
    CGAL/Cartesian/basic_constructions_2.h \
    CGAL/Cartesian/basic_constructions_3.h \
    CGAL/Cartesian/Cartesian_base.h \
    CGAL/Cartesian/Circle_2.h \
    CGAL/Cartesian/Circle_3.h \
    CGAL/Cartesian/ConicCPA2.h \
    CGAL/Cartesian/Data_accessor_2.h \
    CGAL/Cartesian/Direction_2.h \
    CGAL/Cartesian/Direction_3.h \
    CGAL/Cartesian/ft_constructions_2.h \
    CGAL/Cartesian/ft_constructions_3.h \
    CGAL/Cartesian/function_objects.h \
    CGAL/Cartesian/Iso_cuboid_3.h \
    CGAL/Cartesian/Iso_rectangle_2.h \
    CGAL/Cartesian/Line_2.h \
    CGAL/Cartesian/Line_3.h \
    CGAL/Cartesian/line_constructions_2.h \
    CGAL/Cartesian/MatrixC33.h \
    CGAL/Cartesian/Plane_3.h \
    CGAL/Cartesian/plane_constructions_3.h \
    CGAL/Cartesian/Point_2.h \
    CGAL/Cartesian/Point_3.h \
    CGAL/Cartesian/point_constructions_2.h \
    CGAL/Cartesian/point_constructions_3.h \
    CGAL/Cartesian/predicates_on_directions_2.h \
    CGAL/Cartesian/predicates_on_planes_3.h \
    CGAL/Cartesian/predicates_on_points_2.h \
    CGAL/Cartesian/predicates_on_points_3.h \
    CGAL/Cartesian/Ray_2.h \
    CGAL/Cartesian/Ray_3.h \
    CGAL/Cartesian/Rotation_rep_2.h \
    CGAL/Cartesian/Scaling_rep_2.h \
    CGAL/Cartesian/Scaling_rep_3.h \
    CGAL/Cartesian/Segment_2.h \
    CGAL/Cartesian/Segment_3.h \
    CGAL/Cartesian/solve_3.h \
    CGAL/Cartesian/Sphere_3.h \
    CGAL/Cartesian/Tetrahedron_3.h \
    CGAL/Cartesian/Translation_rep_2.h \
    CGAL/Cartesian/Translation_rep_3.h \
    CGAL/Cartesian/Triangle_2.h \
    CGAL/Cartesian/Triangle_3.h \
    CGAL/Cartesian/Vector_2.h \
    CGAL/Cartesian/Vector_3.h \
    CGAL/Cartesian/Weighted_point_2.h \
    CGAL/Cartesian/Weighted_point_3.h \
#    CGAL/Circular_kernel_2/Circular_arc_2.h \
#    CGAL/Circular_kernel_2/Circular_arc_point_2.h \
#    CGAL/Circular_kernel_2/function_objects_on_circle_2.h \
#    CGAL/Circular_kernel_2/function_objects_on_line_2.h \
#    CGAL/Circular_kernel_2/function_objects_polynomial_circular.h \
#    CGAL/Circular_kernel_2/interface_macros.h \
#    CGAL/Circular_kernel_2/internal_functions_on_circle_2.h \
#    CGAL/Circular_kernel_2/internal_functions_on_circular_arc_2.h \
#    CGAL/Circular_kernel_2/internal_functions_on_line_2.h \
#    CGAL/Circular_kernel_2/internal_functions_on_line_arc_2.h \
#    CGAL/Circular_kernel_2/intersection_line_2_circle_2_map.h \
#    CGAL/Circular_kernel_2/Intersection_traits.h \
#    CGAL/Circular_kernel_2/Line_arc_2.h \
#    CGAL/Circular_kernel_3/Circular_arc_3.h \
#    CGAL/Circular_kernel_3/Circular_arc_point_3.h \
#    CGAL/Circular_kernel_3/function_objects_polynomial_sphere.h \
#    CGAL/Circular_kernel_3/get_equation_object_on_curved_kernel_3.h \
#    CGAL/Circular_kernel_3/interface_macros.h \
#    CGAL/Circular_kernel_3/internal_function_compare_spherical_kernel.h \
#    CGAL/Circular_kernel_3/internal_function_compare_to_right_spherical_kernel.h \
#    CGAL/Circular_kernel_3/internal_function_has_on_spherical_kernel.h \
#    CGAL/Circular_kernel_3/internal_functions_on_circle_3.h \
#    CGAL/Circular_kernel_3/internal_functions_on_circular_arc_3.h \
#    CGAL/Circular_kernel_3/internal_functions_on_circular_arc_point_3.h \
#    CGAL/Circular_kernel_3/internal_functions_on_line_3.h \
#    CGAL/Circular_kernel_3/internal_functions_on_line_arc_3.h \
#    CGAL/Circular_kernel_3/internal_functions_on_plane_3.h \
#    CGAL/Circular_kernel_3/internal_functions_on_sphere_3.h \
#    CGAL/Circular_kernel_3/Intersection_traits.h \
#    CGAL/Circular_kernel_3/Line_arc_3.h \
#    CGAL/Circulator/Circulator_adapters.h \
#    CGAL/Circulator/Circulator_concepts.h \
#    CGAL/Circulator/Safe_circulator_from_iterator.h \
    CGAL/Classification/Feature/Cluster_mean_of_feature.h \
    CGAL/Classification/Feature/Cluster_size.h \
    CGAL/Classification/Feature/Cluster_variance_of_feature.h \
    CGAL/Classification/Feature/Cluster_vertical_extent.h \
    CGAL/Classification/Feature/Color_channel.h \
    CGAL/Classification/Feature/Distance_to_plane.h \
    CGAL/Classification/Feature/Echo_scatter.h \
    CGAL/Classification/Feature/Eigen.h \
    CGAL/Classification/Feature/Eigenvalue.h \
    CGAL/Classification/Feature/Elevation.h \
    CGAL/Classification/Feature/Gradient_of_feature.h \
    CGAL/Classification/Feature/Hsv.h \
    CGAL/Classification/Feature/Simple_feature.h \
#    CGAL/Classification/Feature/Vertical_dispersion.h \
#    CGAL/Classification/Feature/Verticality.h \
#    CGAL/Classification/internal/auxiliary/random-forest/common-libraries.hpp \
#    CGAL/Classification/internal/auxiliary/random-forest/forest.hpp \
#    CGAL/Classification/internal/auxiliary/random-forest/node-gini.hpp \
#    CGAL/Classification/internal/auxiliary/random-forest/node.hpp \
#    CGAL/Classification/internal/auxiliary/random-forest/tree.hpp \
#    CGAL/Classification/internal/auxiliary/dataview.h \
#    CGAL/Classification/internal/verbosity.h \
#    CGAL/Classification/classify.h \
#    CGAL/Classification/Cluster.h \
#    CGAL/Classification/Color.h \
#    CGAL/Classification/compressed_float.h \
#    CGAL/Classification/ETHZ_random_forest_classifier.h \
#    CGAL/Classification/Evaluation.h \
#    CGAL/Classification/Feature_base.h \
#    CGAL/Classification/Feature_set.h \
#    CGAL/Classification/Image.h \
#    CGAL/Classification/Label.h \
#    CGAL/Classification/Label_set.h \
#    CGAL/Classification/Local_eigen_analysis.h \
#    CGAL/Classification/Mesh_feature_generator.h \
#    CGAL/Classification/Mesh_neighborhood.h \
#    CGAL/Classification/OpenCV_random_forest_classifier.h \
#    CGAL/Classification/Planimetric_grid.h \
#    CGAL/Classification/Point_set_feature_generator.h \
#    CGAL/Classification/Point_set_neighborhood.h \
#    CGAL/Classification/property_maps.h \
#    CGAL/Classification/Sum_of_weighted_features_classifier.h \
#    CGAL/Cone_spanners_2/Less_by_direction_2.h \
#    CGAL/Cone_spanners_2/Plane_scan_tree.h \
#    CGAL/Cone_spanners_2/Plane_scan_tree_impl.h \
#    CGAL/constructions/constructions_for_voronoi_intersection_cartesian_2_3.h \
#    CGAL/constructions/kernel_ftC2.h \
#    CGAL/constructions/kernel_ftC3.h \
#    CGAL/constructions/Polygon_offset_cons_ftC2.h \
#    CGAL/constructions/Straight_skeleton_cons_ftC2.h \
#    CGAL/Convex_decomposition_3/Edge_sorter.h \
#    CGAL/Convex_decomposition_3/External_structure_builder.h \
#    CGAL/Convex_decomposition_3/Insert_vertex_into_edge.h \
#    CGAL/Convex_decomposition_3/is_reflex_sedge.h \
#    CGAL/Convex_decomposition_3/Ray_hit_generator.h \
#    CGAL/Convex_decomposition_3/Ray_hit_generator2.h \
#    CGAL/Convex_decomposition_3/Reflex_edge_searcher.h \
#    CGAL/Convex_decomposition_3/Reflex_vertex_searcher.h \
#    CGAL/Convex_decomposition_3/SFace_separator.h \
#    CGAL/Convex_decomposition_3/Single_wall_creator.h \
#    CGAL/Convex_decomposition_3/Single_wall_creator2.h \
#    CGAL/Convex_decomposition_3/Single_wall_creator3.h \
#    CGAL/Convex_decomposition_3/SM_walls.h \
#    CGAL/Convex_decomposition_3/YVertical_wall_builder.h \
#    CGAL/Convex_hull_2/ch_akl_toussaint_impl.h \
#    CGAL/Convex_hull_2/ch_assertions.h \
#    CGAL/Convex_hull_2/ch_bykat_impl.h \
#    CGAL/Convex_hull_2/ch_eddy_impl.h \
#    CGAL/Convex_hull_2/ch_graham_andrew_impl.h \
#    CGAL/Convex_hull_2/ch_jarvis_impl.h \
#    CGAL/Convex_hull_2/ch_melkman_impl.h \
#    CGAL/Convex_hull_2/ch_selected_extreme_points_2_impl.h \
#    CGAL/Convex_hull_2/convexity_check_2_impl.h \
#    CGAL/Convex_hull_3/dual/Convex_hull_traits_dual_2.h \
#    CGAL/Convex_hull_3/dual/Convex_hull_traits_dual_3.h \
#    CGAL/Convex_hull_3/dual/halfspace_intersection_3.h \
#    CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h \
#    CGAL/Convex_hull_3/dual/interior_polyhedron_3.h \
#    CGAL/Convex_hull_3/dual/predicates.h \
#    CGAL/CORE/poly/Curves.h \
#    CGAL/CORE/poly/Poly.h \
#    CGAL/CORE/poly/Sturm.h \
#    CGAL/CORE/BigFloat.h \
#    CGAL/CORE/BigFloat_impl.h \
#    CGAL/CORE/BigFloatRep.h \
#    CGAL/CORE/BigInt.h \
#    CGAL/CORE/BigRat.h \
#    CGAL/CORE/Config.h \
#    CGAL/CORE/CORE.h \
#    CGAL/CORE/CoreAux.h \
#    CGAL/CORE/CoreAux_impl.h \
#    CGAL/CORE/CoreDefs.h \
#    CGAL/CORE/CoreDefs_impl.h \
#    CGAL/CORE/CoreIO_impl.h \
#    CGAL/CORE/Expr.h \
#    CGAL/CORE/Expr_impl.h \
#    CGAL/CORE/ExprRep.h \
#    CGAL/CORE/extLong.h \
#    CGAL/CORE/extLong_impl.h \
#    CGAL/CORE/Filter.h \
#    CGAL/CORE/Gmp.h \
#    CGAL/CORE/Gmp_impl.h \
#    CGAL/CORE/Impl.h \
#    CGAL/CORE/linearAlgebra.h \
#    CGAL/CORE/MemoryPool.h \
#    CGAL/CORE/Promote.h \
#    CGAL/CORE/Real.h \
#    CGAL/CORE/Real_impl.h \
#    CGAL/CORE/RealRep.h \
#    CGAL/CORE/RefCount.h \
#    CGAL/CORE/Timer.h \
#    CGAL/Curved_kernel_via_analysis_2/gfx/Curve_renderer_2.h \
#    CGAL/Curved_kernel_via_analysis_2/gfx/Curve_renderer_internals.h \
#    CGAL/Curved_kernel_via_analysis_2/gfx/Curve_renderer_traits.h \
#    CGAL/Curved_kernel_via_analysis_2/gfx/Subdivision_1.h \
#    CGAL/Curved_kernel_via_analysis_2/gfx/Subdivision_2.h \
#    CGAL/Curved_kernel_via_analysis_2/test/simple_models.h \
#    CGAL/Curved_kernel_via_analysis_2/Arc_2.h \
#    CGAL/Curved_kernel_via_analysis_2/Curve_interval_arcno_cache.h \
#    CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h \
#    CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h \
#    CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h \
#    CGAL/Curved_kernel_via_analysis_2/Fig_stream_Curve_renderer_2.h \
#    CGAL/Curved_kernel_via_analysis_2/Filtered_curved_kernel_via_analysis_2_impl.h \
#    CGAL/Curved_kernel_via_analysis_2/Generic_arc_2.h \
#    CGAL/Curved_kernel_via_analysis_2/Generic_point_2.h \
#    CGAL/Curved_kernel_via_analysis_2/Make_x_monotone_2.h \
#    CGAL/Curved_kernel_via_analysis_2/Non_x_monotone_arc_2.h \
#    CGAL/Curved_kernel_via_analysis_2/Point_2.h \
#    CGAL/Curved_kernel_via_analysis_2/Sweep_curves_adapter_2.h \
#    CGAL/Envelope_2/Env_divide_and_conquer_2.h \
#    CGAL/Envelope_2/Env_divide_and_conquer_2_impl.h \
#    CGAL/Envelope_3/Env_plane_traits_3_functions.h \
#    CGAL/Envelope_3/Envelope_base.h \
#    CGAL/Envelope_3/Envelope_diagram_on_surface_2.h \
#    CGAL/Envelope_3/Envelope_divide_and_conquer_3.h \
#    CGAL/Envelope_3/Envelope_element_visitor_3.h \
#    CGAL/Envelope_3/Envelope_overlay_2.h \
#    CGAL/Envelope_3/Envelope_overlay_functor.h \
#    CGAL/Envelope_3/Envelope_pm_dcel.h \
#    CGAL/Envelope_3/set_dividors.h \
#    CGAL/export/CGAL.h \
#    CGAL/export/CORE.h \
#    CGAL/export/helpers.h \
#    CGAL/export/ImageIO.h \
#    CGAL/export/Qt.h \
#    CGAL/Filtered_bbox_circular_kernel_2/bbox_filtered_predicates.h \
#    CGAL/Filtered_bbox_circular_kernel_2/interface_macros.h \
#    CGAL/Filtered_kernel/Cartesian_coordinate_iterator_2.h \
#    CGAL/Filtered_kernel/Cartesian_coordinate_iterator_3.h \
#    CGAL/GMP/Gmpfi_type.h \
#    CGAL/GMP/Gmpfi_type_static.h \
#    CGAL/GMP/Gmpfr_type.h \
#    CGAL/GMP/Gmpfr_type_static.h \
#    CGAL/GMP/Gmpq_type.h \
#    CGAL/GMP/Gmpz_type.h \
#    CGAL/GMP/Gmpzf_type.h \
#    CGAL/Homogeneous/Aff_transformationH2.h \
#    CGAL/Homogeneous/Aff_transformationH3.h \
#    CGAL/Homogeneous/basic_constructionsH2.h \
#    CGAL/Homogeneous/basic_constructionsH3.h \
#    CGAL/Homogeneous/CircleH2.h \
#    CGAL/Homogeneous/ConicHPA2.h \
#    CGAL/Homogeneous/Data_accessorH2.h \
#    CGAL/Homogeneous/DirectionH2.h \
#    CGAL/Homogeneous/DirectionH3.h \
#    CGAL/Homogeneous/distance_predicatesH2.h \
#    CGAL/Homogeneous/distance_predicatesH3.h \
#    CGAL/Homogeneous/function_objects.h \
#    CGAL/Homogeneous/Homogeneous_base.h \
#    CGAL/Homogeneous/Iso_cuboidH3.h \
#    CGAL/Homogeneous/Iso_rectangleH2.h \
#    CGAL/Homogeneous/LineH2.h \
#    CGAL/Homogeneous/PlaneH3.h \
#    CGAL/Homogeneous/PointH2.h \
#    CGAL/Homogeneous/PointH3.h \
#    CGAL/Homogeneous/predicates_on_directionsH2.h \
#    CGAL/Homogeneous/predicates_on_pointsH2.h \
#    CGAL/Homogeneous/predicates_on_pointsH3.h \
#    CGAL/Homogeneous/RayH3.h \
#    CGAL/Homogeneous/SphereH3.h \
#    CGAL/Homogeneous/VectorH2.h \
#    CGAL/Homogeneous/VectorH3.h \
#    CGAL/Homogeneous/Weighted_point_2.h \
#    CGAL/Homogeneous/Weighted_point_3.h \
    CGAL/ImageIO/analyze.h \
    CGAL/ImageIO/analyze_impl.h \
    CGAL/ImageIO/bmp.h \
    CGAL/ImageIO/bmp_impl.h \
    CGAL/ImageIO/bmpendian.h \
    CGAL/ImageIO/bmpendian_impl.h \
    CGAL/ImageIO/bmpread.h \
    CGAL/ImageIO/bmpread_impl.h \
    CGAL/ImageIO/bmptypes.h \
    CGAL/ImageIO/convert.h \
    CGAL/ImageIO/convert_impl.h \
    CGAL/ImageIO/fgetns.h \
    CGAL/ImageIO/fgetns_impl.h \
    CGAL/ImageIO/gif.h \
    CGAL/ImageIO/gif_impl.h \
    CGAL/ImageIO/gis.h \
    CGAL/ImageIO/gis_impl.h \
    CGAL/ImageIO/inr.h \
    CGAL/ImageIO/inr_impl.h \
    CGAL/ImageIO/iris.h \
    CGAL/ImageIO/iris_impl.h \
    CGAL/ImageIO/mincio.h \
    CGAL/ImageIO/mincio_impl.h \
    CGAL/ImageIO/pnm.h \
    CGAL/ImageIO/pnm_impl.h \
    CGAL/ImageIO/recbuffer.h \
    CGAL/ImageIO/recbuffer_impl.h \
    CGAL/ImageIO/recline.h \
    CGAL/ImageIO/recline_impl.h \
    CGAL/ImageIO/reech4x4.h \
    CGAL/ImageIO/reech4x4_impl.h \
    CGAL/ImageIO/typedefs.h \
#    CGAL/internal/AABB_tree/AABB_drawing_traits.h \
#    CGAL/internal/AABB_tree/AABB_node.h \
#    CGAL/internal/AABB_tree/AABB_ray_intersection.h \
#    CGAL/internal/AABB_tree/AABB_search_tree.h \
#    CGAL/internal/AABB_tree/AABB_traversal_traits.h \
#    CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h \
#    CGAL/internal/AABB_tree/Is_ray_intersection_geomtraits.h \
#    CGAL/internal/AABB_tree/Primitive_helper.h \
#    CGAL/internal/AFSR/construct_polyhedron.h \
#    CGAL/internal/AFSR/construct_surface_2.h \
#    CGAL/internal/AFSR/orient.h \
#    CGAL/internal/AFSR/Surface_face_base_2.h \
#    CGAL/internal/AFSR/Surface_vertex_base_2.h \
#    CGAL/internal/AFSR/write_triple_indices.h \
#    CGAL/internal/auxiliary/graph.h \
#    CGAL/internal/boost/array_binary_tree.hpp \
#    CGAL/internal/boost/function_property_map.hpp \
#    CGAL/internal/boost/mutable_heap.hpp \
#    CGAL/internal/boost/mutable_queue.hpp \
#    CGAL/internal/corefinement/Combinatorial_map_for_corefinement.h \
#    CGAL/internal/corefinement/Combinatorial_map_output_builder.h \
#    CGAL/internal/corefinement/connected_components.h \
#    CGAL/internal/corefinement/intersection_coplanar_triangles_3.h \
#    CGAL/internal/corefinement/intersection_triangle_segment_3.h \
#    CGAL/internal/corefinement/intersection_triangle_segment_3_coplanar.h \
#    CGAL/internal/corefinement/Polyhedra_output_builder.h \
#    CGAL/internal/corefinement/Polyhedron_constness_types.h \
#    CGAL/internal/corefinement/predicates.h \
#    CGAL/internal/corefinement/utils.h \
#    CGAL/internal/Intersections_3/Bbox_3_Bbox_3_do_intersect.h \
#    CGAL/internal/Intersections_3/Bbox_3_Line_3_do_intersect.h \
#    CGAL/internal/Intersections_3/Bbox_3_Plane_3_do_intersect.h \
#    CGAL/internal/Intersections_3/Bbox_3_Ray_3_do_intersect.h \
#    CGAL/internal/Intersections_3/Bbox_3_Segment_3_do_intersect.h \
#    CGAL/internal/Intersections_3/Bbox_3_Sphere_3_do_intersect.h \
#    CGAL/internal/Intersections_3/Bbox_3_Triangle_3_do_intersect.h \
#    CGAL/internal/Intersections_3/Triangle_3_Line_3_intersection.h \
#    CGAL/internal/Intersections_3/Triangle_3_Ray_3_intersection.h \
#    CGAL/internal/Intersections_3/Triangle_3_Segment_3_intersection.h \
#    CGAL/internal/Mesh_3/Boundary_of_subdomain_of_complex_3_in_triangulation_3_to_off.h \
#    CGAL/internal/Mesh_3/check_weights.h \
#    CGAL/internal/Mesh_3/get_index.h \
#    CGAL/internal/Mesh_3/Graph_manipulations.h \
#    CGAL/internal/Mesh_3/Handle_IO_for_pair_of_int.h \
#    CGAL/internal/Mesh_3/helpers.h \
#    CGAL/internal/Mesh_3/indices_management.h \
#    CGAL/internal/Static_filters/Angle_3.h \
#    CGAL/internal/Static_filters/Collinear_3.h \
#    CGAL/internal/Static_filters/Compare_squared_radius_3.h \
#    CGAL/internal/Static_filters/Compare_weighted_squared_radius_3.h \
#    CGAL/internal/Static_filters/Compare_x_2.h \
#    CGAL/internal/Static_filters/Compare_y_2.h \
#    CGAL/internal/Static_filters/Compare_y_at_x_2.h \
#    CGAL/internal/Static_filters/Coplanar_orientation_3.h \
#    CGAL/internal/Static_filters/Coplanar_side_of_bounded_circle_3.h \
#    CGAL/internal/Static_filters/Do_intersect_2.h \
#    CGAL/internal/Static_filters/Do_intersect_3.h \
#    CGAL/internal/Static_filters/Equal_2.h \
#    CGAL/internal/Static_filters/Equal_3.h \
#    CGAL/internal/Static_filters/Is_degenerate_3.h \
#    CGAL/internal/Static_filters/Orientation_2.h \
#    CGAL/internal/Static_filters/Orientation_3.h \
#    CGAL/internal/Static_filters/Periodic_2_orientation_2.h \
#    CGAL/internal/Static_filters/Periodic_2_side_of_oriented_circle_2.h \
#    CGAL/internal/Static_filters/Periodic_3_orientation_3.h \
#    CGAL/internal/Static_filters/Periodic_3_power_side_of_oriented_power_sphere_3.h \
#    CGAL/internal/Static_filters/Periodic_3_side_of_oriented_sphere_3.h \
#    CGAL/internal/Static_filters/Power_side_of_oriented_power_sphere_3.h \
#    CGAL/internal/Static_filters/Side_of_oriented_circle_2.h \
#    CGAL/internal/Static_filters/Side_of_oriented_sphere_3.h \
#    CGAL/internal/Static_filters/Static_filter_error.h \
#    CGAL/internal/Static_filters/Static_filters.h \
#    CGAL/internal/Static_filters/tools.h \
#    CGAL/internal/Surface_mesh_deformation/Spokes_and_rims_iterator.h \
#    CGAL/internal/Surface_mesh_segmentation/AABB_traits.h \
#    CGAL/internal/Surface_mesh_segmentation/AABB_traversal_traits.h \
#    CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h \
#    CGAL/internal/Surface_mesh_segmentation/Disk_samplers.h \
#    CGAL/internal/Surface_mesh_segmentation/Expectation_maximization.h \
#    CGAL/internal/Surface_mesh_segmentation/Filters.h \
#    CGAL/internal/Surface_mesh_segmentation/K_means_clustering.h \
#    CGAL/internal/Surface_mesh_segmentation/SDF_calculation.h \
#    CGAL/internal/Surface_mesh_segmentation/Surface_mesh_segmentation.h \
#    CGAL/internal/Surface_mesh_skeletonization/Curve_skeleton.h \
#    CGAL/internal/Surface_mesh_skeletonization/Debug.h \
#    CGAL/internal/Surface_mesh_skeletonization/Detect_degeneracy.h \
#    CGAL/internal/TDS_2/Edge_hash_function.h \
#    CGAL/internal/TDS_2/edge_list.h \
#    CGAL/internal/Triangulation/Dummy_TDS.h \
#    CGAL/internal/Triangulation/Triangulation_ds_iterators.h \
#    CGAL/internal/Triangulation/utilities.h \
#    CGAL/internal/bounded_priority_queue.h \
#    CGAL/internal/canonicalize_helper.h \
#    CGAL/internal/Classification_type.h \
#    CGAL/internal/Combination_enumerator.h \
#    CGAL/internal/Combinatorial_map_copy_functors.h \
#    CGAL/internal/Combinatorial_map_group_functors.h \
#    CGAL/internal/Combinatorial_map_internal_functors.h \
#    CGAL/internal/Combinatorial_map_sewable.h \
#    CGAL/internal/Combinatorial_map_utility.h \
#    CGAL/internal/Combinatorial_map_utility_novariadic.h \
#    CGAL/internal/container_fwd_fixed.hpp \
#    CGAL/internal/Delaunay_triangulation_hierarchy_3.h \
#    CGAL/internal/deprecation_warning.h \
#    CGAL/internal/disable_deprecation_warnings_and_errors.h \
#    CGAL/internal/Dummy_tds_3.h \
#    CGAL/internal/enable_third_party_libraries.h \
#    CGAL/internal/Exact_type_selector.h \
#    CGAL/internal/Functor_with_offset_points_adaptor_2.h \
#    CGAL/internal/Functor_with_offset_points_adaptor_3.h \
#    CGAL/internal/Functor_with_offset_weighted_points_adaptor_3.h \
#    CGAL/internal/Generalized_map_group_functors.h \
#    CGAL/internal/Generalized_map_internal_functors.h \
#    CGAL/internal/Generalized_map_sewable.h \
#    CGAL/internal/Generic_random_point_generator.h \
#    CGAL/internal/Get_dimension_tag.h \
#    CGAL/internal/Has_boolean_tags.h \
#    CGAL/internal/Has_member_visited.h \
#    CGAL/internal/Has_nested_type_Bare_point.h \
#    CGAL/internal/info_check.h \
#    CGAL/internal/K_neighbor_search.h \
#    CGAL/internal/Lazy_alpha_nt_2.h \
#    CGAL/internal/Lazy_alpha_nt_3.h \
#    CGAL/internal/Parallel_callback.h \
#    CGAL/internal/Periodic_2_construct_point_2.h \
#    CGAL/internal/Periodic_2_Delaunay_triangulation_filtered_traits_2.h \
#    CGAL/internal/Periodic_2_Delaunay_triangulation_statically_filtered_traits_2.h \
#    CGAL/internal/Periodic_2_triangulation_filtered_traits_2.h \
#    CGAL/internal/Periodic_2_triangulation_statically_filtered_traits_2.h \
#    CGAL/internal/Periodic_3_construct_point_3.h \
#    CGAL/internal/Periodic_3_construct_weighted_point_3.h \
#    CGAL/internal/Periodic_3_Delaunay_triangulation_filtered_traits_3.h \
#    CGAL/internal/Periodic_3_Delaunay_triangulation_remove_traits_3.h \
#    CGAL/internal/Periodic_3_Delaunay_triangulation_statically_filtered_traits_3.h \
#    CGAL/internal/Periodic_3_regular_triangulation_dummy_288.h \
#    CGAL/internal/Periodic_3_regular_triangulation_filtered_traits_3.h \
#    CGAL/internal/Periodic_3_regular_triangulation_remove_traits_3.h \
#    CGAL/internal/Periodic_3_regular_triangulation_statically_filtered_traits_3.h \
#    CGAL/internal/Periodic_3_triangulation_dummy_36.h \
#    CGAL/internal/Periodic_3_triangulation_filtered_traits_3.h \
#    CGAL/internal/Periodic_3_triangulation_iterators_3.h \
#    CGAL/internal/Periodic_3_triangulation_remove_traits_3.h \
#    CGAL/internal/Periodic_3_triangulation_statically_filtered_traits_3.h \
#    CGAL/internal/Projection_traits_3.h \
#    CGAL/internal/Robust_periodic_weighted_circumcenter_traits_3.h \
#    CGAL/internal/Static_or_dynamic_array.h \
#    CGAL/internal/Transform_coordinates_traits_3.h \
#    CGAL/internal/Triangulation_2_filtered_projection_traits_3.h \
#    CGAL/internal/Triangulation_2_projection_traits_base_3.h \
#    CGAL/internal/Triangulation_ds_circulators_3.h \
#    CGAL/internal/Triangulation_ds_iterators_3.h \
#    CGAL/Interpolation/internal/helpers.h \
#    CGAL/Intersections_2/Triangle_2_Triangle_2_intersection_impl.h \
#    CGAL/Intersections_3/intersection_3_1_impl.h \
    CGAL/IO/Alpha_shape_3_VRML_2_ostream.h \
    CGAL/IO/alpha_shape_geomview_ostream_3.h \
    CGAL/IO/Arr_iostream.h \
    CGAL/IO/Arr_text_formatter.h \
    CGAL/IO/Arr_with_history_2_reader.h \
    CGAL/IO/Arr_with_history_2_writer.h \
    CGAL/IO/Arr_with_history_iostream.h \
    CGAL/IO/Arr_with_history_text_formatter.h \
    CGAL/IO/Arrangement_2_reader.h \
    CGAL/IO/Arrangement_2_writer.h \
    CGAL/IO/binary_file_io.h \
    CGAL/IO/Color.h \
    CGAL/IO/Color_impl.h \
    CGAL/IO/Complex_2_in_triangulation_3_file_writer.h \
    CGAL/IO/Complex_2_in_triangulation_3_polyhedron_builder.h \
    CGAL/IO/Complex_2_in_triangulation_3_to_medit.h \
    CGAL/IO/Complex_2_in_triangulation_3_to_vtk.h \
    CGAL/IO/Complex_3_in_triangulation_3_to_vtk.h \
    CGAL/IO/Dxf_bsop_reader.h \
    CGAL/IO/Dxf_reader.h \
    CGAL/IO/Dxf_reader_doubles.h \
    CGAL/IO/Dxf_stream.h \
    CGAL/IO/Dxf_variant_reader.h \
    CGAL/IO/Dxf_writer.h \
    CGAL/IO/facets_in_complex_2_to_triangle_mesh.h \
    CGAL/IO/facets_in_complex_3_to_triangle_mesh.h \
    CGAL/IO/Fig_stream.h \
    CGAL/IO/Fig_stream_Conic_arc_2.h \
    CGAL/IO/File_avizo.h \
    CGAL/IO/File_binary_mesh_3.h \
    CGAL/IO/File_header_extended_OFF.h \
    CGAL/IO/File_header_extended_OFF_impl.h \
    CGAL/IO/File_header_OFF.h \
    CGAL/IO/File_header_OFF_impl.h \
    CGAL/IO/File_maya.h \
    CGAL/IO/File_medit.h \
    CGAL/IO/File_poly.h \
    CGAL/IO/File_scanner_OFF.h \
    CGAL/IO/File_scanner_OFF_impl.h \
    CGAL/IO/File_tetgen.h \
    CGAL/IO/File_writer_inventor.h \
    CGAL/IO/File_writer_inventor_impl.h \
    CGAL/IO/File_writer_OFF.h \
    CGAL/IO/File_writer_OFF_impl.h \
    CGAL/IO/File_writer_VRML_2.h \
    CGAL/IO/File_writer_VRML_2_impl.h \
    CGAL/IO/File_writer_wavefront.h \
    CGAL/IO/File_writer_wavefront_impl.h \
    CGAL/IO/generic_copy_OFF.h \
    CGAL/IO/generic_print_polyhedron.h \
    CGAL/IO/Generic_writer.h \
    CGAL/IO/Geomview_stream.h \
    CGAL/IO/Geomview_stream_impl.h \
    CGAL/IO/Gps_iostream.h \
    CGAL/IO/Inventor_ostream.h \
    CGAL/IO/io.h \
    CGAL/IO/io_impl.h \
    CGAL/IO/io_tags.h \
    CGAL/IO/Istream_iterator.h \
    CGAL/IO/Nef_polyhedron_2_PS_stream.h \
    CGAL/IO/Nef_polyhedron_iostream_3.h \
    CGAL/IO/OBJ_reader.h \
    CGAL/IO/OFF_reader.h \
    CGAL/IO/Ostream_iterator.h \
    CGAL/IO/output_surface_facets_to_polyhedron.h \
    CGAL/IO/output_surface_facets_to_triangle_soup.h \
    CGAL/IO/PLY_reader.h \
    CGAL/IO/PLY_writer.h \
    CGAL/IO/Polyhedron_builder_from_STL.h \
    CGAL/IO/Polyhedron_geomview_ostream.h \
    CGAL/IO/Polyhedron_inventor_ostream.h \
    CGAL/IO/Polyhedron_iostream.h \
    CGAL/IO/Polyhedron_scan_OFF.h \
    CGAL/IO/Polyhedron_VRML_1_ostream.h \
    CGAL/IO/Polyhedron_VRML_2_ostream.h \
    CGAL/IO/print_inventor.h \
    CGAL/IO/print_OFF.h \
    CGAL/IO/print_VRML_1.h \
    CGAL/IO/print_VRML_2.h \
    CGAL/IO/print_wavefront.h \
    CGAL/IO/read_las_points.h \
    CGAL/IO/read_off_points.h \
    CGAL/IO/read_ply_points.h \
    CGAL/IO/read_xyz_points.h \
    CGAL/IO/scan_OFF.h \
    CGAL/IO/Scanner_OFF.h \
    CGAL/IO/STL_reader.h \
    CGAL/IO/STL_writer.h \
    CGAL/IO/Tee_for_output_iterator.h \
    CGAL/IO/Triangulation_geomview_ostream_2.h \
    CGAL/IO/Triangulation_geomview_ostream_3.h \
    CGAL/IO/Triangulation_off_ostream.h \
    CGAL/IO/Triangulation_off_ostream_2.h \
    CGAL/IO/Triangulation_off_ostream_3.h \
    CGAL/IO/Triangulation_ps_stream.h \
    CGAL/IO/Verbose_ostream.h \
    CGAL/IO/VRML_1_ostream.h \
    CGAL/IO/VRML_2_ostream.h \
    CGAL/IO/write_las_points.h \
    CGAL/IO/write_off_points.h \
    CGAL/IO/write_ply_points.h \
    CGAL/IO/write_xyz_points.h \
    CGAL/IO/Writer_OFF.h \
#    CGAL/Kernel/Conic_misc.h \
#    CGAL/Kernel/Dimension_utils.h \
#    CGAL/Kernel/function_objects.h \
#    CGAL/Kernel/global_functions.h \
#    CGAL/Kernel/global_functions_2.h \
#    CGAL/Kernel/global_functions_3.h \
#    CGAL/Kernel/global_functions_internal_2.h \
#    CGAL/Kernel/global_functions_internal_3.h \
#    CGAL/Kernel/interface_macros.h \
#    CGAL/Kernel/mpl.h \
#    CGAL/Kernel/Return_base_tag.h \
#    CGAL/Kernel/Same_uncertainty.h \
#    CGAL/Kernel/solve.h \
#    CGAL/Kernel/Type_equality_wrapper.h \
#    CGAL/Kernel/Type_mapper.h \
#    CGAL/Kernel/Wutils.h \
#    CGAL/Kernel_d/Aff_transformation_d.h \
#    CGAL/Kernel_d/Aff_transformationCd.h \
#    CGAL/Kernel_d/Aff_transformationHd.h \
#    CGAL/Kernel_d/Cartesian_const_iterator_d.h \
#    CGAL/Kernel_d/Cartesian_converter_d.h \
#    CGAL/Kernel_d/debug.h \
#    CGAL/Kernel_d/Direction_d.h \
#    CGAL/Kernel_d/DirectionCd.h \
#    CGAL/Kernel_d/DirectionCd_impl.h \
#    CGAL/Kernel_d/DirectionHd.h \
#    CGAL/Kernel_d/DirectionHd_impl.h \
#    CGAL/Kernel_d/function_objects.h \
#    CGAL/Kernel_d/function_objectsCd.h \
#    CGAL/Kernel_d/function_objectsHd.h \
#    CGAL/Kernel_d/Hyperplane_d.h \
#    CGAL/Kernel_d/HyperplaneCd.h \
#    CGAL/Kernel_d/HyperplaneCd_impl.h \
#    CGAL/Kernel_d/HyperplaneHd.h \
#    CGAL/Kernel_d/HyperplaneHd_impl.h \
#    CGAL/Kernel_d/Interface_classes.h \
#    CGAL/Kernel_d/interface_macros_d.h \
#    CGAL/Kernel_d/intersection_objects_d.h \
#    CGAL/Kernel_d/intersection_objectsCd.h \
#    CGAL/Kernel_d/intersection_objectsHd.h \
#    CGAL/Kernel_d/Interval_linear_algebra.h \
#    CGAL/Kernel_d/Iso_box_d.h \
#    CGAL/Kernel_d/Kernel_classesCd.h \
#    CGAL/Kernel_d/Kernel_classesHd.h \
#    CGAL/Kernel_d/Line_d.h \
#    CGAL/Kernel_d/Line_d_impl.h \
#    CGAL/Kernel_d/Linear_algebraCd_impl.h \
#    CGAL/Kernel_d/Linear_algebraHd_impl.h \
#    CGAL/Kernel_d/Matrix__.h \
#    CGAL/Kernel_d/Pair_d.h \
#    CGAL/Kernel_d/Point_d.h \
#    CGAL/Kernel_d/PointCd.h \
#    CGAL/Kernel_d/PointCd_impl.h \
#    CGAL/Kernel_d/PointHd.h \
#    CGAL/Kernel_d/PointHd_impl.h \
#    CGAL/Kernel_d/Ray_d.h \
#    CGAL/Kernel_d/Segment_d.h \
#    CGAL/Kernel_d/simple_objects.h \
#    CGAL/Kernel_d/Sphere_d.h \
#    CGAL/Kernel_d/Tuple_d.h \
#    CGAL/Kernel_d/Vector__.h \
#    CGAL/Kernel_d/Vector_d.h \
#    CGAL/Kernel_d/VectorCd.h \
#    CGAL/Kernel_d/VectorCd_impl.h \
#    CGAL/Kernel_d/VectorHd.h \
#    CGAL/Kernel_d/VectorHd_impl.h \
#    CGAL/license/Polygon_mesh_processing/Compute_normal.h \
#    CGAL/license/Polygon_mesh_processing/connected_components.h \
#    CGAL/license/Polygon_mesh_processing/core.h \
#    CGAL/license/Polygon_mesh_processing/corefinement.h \
#    CGAL/license/Polygon_mesh_processing/detect_features.h \
#    CGAL/license/Polygon_mesh_processing/distance.h \
#    CGAL/license/Polygon_mesh_processing/measure.h \
#    CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h \
#    CGAL/license/Polygon_mesh_processing/miscellaneous.h \
#    CGAL/license/Polygon_mesh_processing/orientation.h \
#    CGAL/license/Polygon_mesh_processing/predicate.h \
#    CGAL/license/Polygon_mesh_processing/repair.h \
#    CGAL/license/AABB_tree.h \
#    CGAL/license/Advancing_front_surface_reconstruction.h \
#    CGAL/license/Alpha_shapes_2.h \
#    CGAL/license/Alpha_shapes_3.h \
#    CGAL/license/Apollonius_graph_2.h \
#    CGAL/license/Arrangement_on_surface_2.h \
#    CGAL/license/Barycentric_coordinates_2.h \
#    CGAL/license/Boolean_set_operations_2.h \
#    CGAL/license/Bounding_volumes.h \
#    CGAL/license/Box_intersection_d.h \
#    CGAL/license/Circular_kernel_2.h \
#    CGAL/license/Circular_kernel_3.h \
#    CGAL/license/Classification.h \
#    CGAL/license/Cone_spanners_2.h \
#    CGAL/license/Convex_decomposition_3.h \
#    CGAL/license/Convex_hull_2.h \
#    CGAL/license/Convex_hull_3.h \
#    CGAL/license/Convex_hull_d.h \
#    CGAL/license/Envelope_2.h \
#    CGAL/license/Envelope_3.h \
#    CGAL/license/GraphicsView.h \
#    CGAL/license/Inscribed_areas.h \
#    CGAL/license/Interpolation.h \
#    CGAL/license/Interval_skip_list.h \
#    CGAL/license/Jet_fitting_3.h \
#    CGAL/license/lgpl.h \
#    CGAL/license/Matrix_search.h \
#    CGAL/license/Mesh_2.h \
#    CGAL/license/Mesh_3.h \
#    CGAL/license/Minkowski_sum_2.h \
#    CGAL/license/Minkowski_sum_3.h \
#    CGAL/license/Nef_2.h \
#    CGAL/license/Nef_3.h \
#    CGAL/license/Nef_S2.h \
#    CGAL/license/Optimal_transportation_reconstruction_2.h \
#    CGAL/license/Partition_2.h \
#    CGAL/license/Periodic_2_triangulation_2.h \
#    CGAL/license/Periodic_3_mesh_3.h \
#    CGAL/license/Periodic_3_triangulation_3.h \
#    CGAL/license/Point_set_2.h \
#    CGAL/license/Point_set_3.h \
#    CGAL/license/Point_set_processing_3.h \
#    CGAL/license/Point_set_shape_detection_3.h \
#    CGAL/license/Poisson_surface_reconstruction_3.h \
#    CGAL/license/Polygon_mesh_processing.h \
#    CGAL/license/Polyhedron.h \
#    CGAL/license/Polyline_simplification_2.h \
#    CGAL/license/Polytope_distance_d.h \
#    CGAL/license/Principal_component_analysis.h \
#    CGAL/license/QP_solver.h \
#    CGAL/license/Ridges_3.h \
#    CGAL/license/Scale_space_reconstruction_3.h \
#    CGAL/license/SearchStructures.h \
#    CGAL/license/Segment_Delaunay_graph_2.h \
#    CGAL/license/Segment_Delaunay_graph_Linf_2.h \
#    CGAL/license/Set_movable_separability_2.h \
#    CGAL/license/Skin_surface_3.h \
#    CGAL/license/Snap_rounding_2.h \
#    CGAL/license/Spatial_searching.h \
#    CGAL/license/Straight_skeleton_2.h \
#    CGAL/license/Stream_lines_2.h \
#    CGAL/license/Surface_mesh.h \
#    CGAL/license/Surface_mesh_deformation.h \
#    CGAL/license/Surface_mesh_parameterization.h \
#    CGAL/license/Surface_mesh_segmentation.h \
#    CGAL/license/Surface_mesh_shortest_path.h \
#    CGAL/license/Surface_mesh_simplification.h \
#    CGAL/license/Surface_mesh_skeletonization.h \
#    CGAL/license/Surface_mesher.h \
#    CGAL/license/Surface_sweep_2.h \
#    CGAL/license/TDS_2.h \
#    CGAL/license/TDS_3.h \
#    CGAL/license/Three.h \
#    CGAL/license/Triangulation.h \
#    CGAL/license/Triangulation_2.h \
#    CGAL/license/Triangulation_3.h \
#    CGAL/license/Visibility_2.h \
#    CGAL/license/Voronoi_diagram_2.h \
    CGAL/Mesh_2/Clusters.h \
    CGAL/Mesh_2/Do_not_refine_edges.h \
    CGAL/Mesh_2/Face_badness.h \
    CGAL/Mesh_2/Lipschitz_sizing_field_2.h \
    CGAL/Mesh_2/Lloyd_move_2.h \
    CGAL/Mesh_2/Mesh_global_optimizer_2.h \
    CGAL/Mesh_2/Mesh_sizing_field.h \
    CGAL/Mesh_2/Output_stream.h \
    CGAL/Mesh_2/Refine_edges.h \
    CGAL/Mesh_2/Refine_edges_visitor.h \
    CGAL/Mesh_2/Refine_edges_with_clusters.h \
    CGAL/Mesh_2/Refine_faces.h \
    CGAL/Mesh_2/Sizing_field_2.h \
    CGAL/Mesh_2/Uniform_sizing_field_2.h \
    CGAL/Mesh_3/experimental/AABB_filtered_projection_traits.h \
    CGAL/Mesh_3/experimental/Facet_topological_criterion_with_adjacency.h \
    CGAL/Mesh_3/experimental/Get_curve_index.h \
    CGAL/Mesh_3/experimental/Get_facet_patch_id.h \
    CGAL/Mesh_3/experimental/Lipschitz_sizing_experimental.h \
    CGAL/Mesh_3/experimental/Lipschitz_sizing_parameters.h \
    CGAL/Mesh_3/experimental/Lipschitz_sizing_polyhedron.h \
    CGAL/Mesh_3/experimental/Sizing_field_minimum.h \
    CGAL/Mesh_3/experimental/Sizing_field_with_aabb_tree.h \
    CGAL/Mesh_3/C3T3_helpers.h \
    CGAL/Mesh_3/Cell_criteria_visitor_with_balls.h \
    CGAL/Mesh_3/Concurrent_mesher_config.h \
    CGAL/Mesh_3/config.h \
    CGAL/Mesh_3/Detect_polylines_in_polyhedra.h \
    CGAL/Mesh_3/Detect_polylines_in_polyhedra_fwd.h \
    CGAL/Mesh_3/dihedral_angle_3.h \
    CGAL/Mesh_3/Dump_c3t3.h \
    CGAL/Mesh_3/Facet_criteria_visitor_with_balls.h \
    CGAL/Mesh_3/Facet_on_same_surface_criterion.h \
    CGAL/Mesh_3/Has_features.h \
    CGAL/Mesh_3/Image_to_labeled_function_wrapper.h \
    CGAL/Mesh_3/Implicit_surface_mesher_visitor.h \
    CGAL/Mesh_3/initialize_triangulation_from_labeled_image.h \
    CGAL/Mesh_3/io_signature.h \
    CGAL/Mesh_3/Lloyd_move.h \
    CGAL/Mesh_3/Mesh_complex_3_in_triangulation_3_base.h \
    CGAL/Mesh_3/Mesh_global_optimizer.h \
    CGAL/Mesh_3/Mesh_sizing_field.h \
    CGAL/Mesh_3/mesh_standard_cell_criteria.h \
    CGAL/Mesh_3/mesh_standard_criteria.h \
    CGAL/Mesh_3/mesh_standard_facet_criteria.h \
    CGAL/Mesh_3/Mesh_surface_cell_base_3.h \
    CGAL/Mesh_3/Mesher_3.h \
    CGAL/Mesh_3/Mesher_level.h \
    CGAL/Mesh_3/Mesher_level_default_implementations.h \
    CGAL/Mesh_3/min_dihedral_angle.h \
    CGAL/Mesh_3/Null_exuder_visitor.h \
    CGAL/Mesh_3/Null_global_optimizer_visitor.h \
    CGAL/Mesh_3/Null_perturber_visitor.h \
    CGAL/Mesh_3/Odt_move.h \
    CGAL/Mesh_3/parameters_defaults.h \
    CGAL/Mesh_3/Poisson_refine_cells_3.h \
    CGAL/Mesh_3/polyhedral_to_labeled_function_wrapper.h \
    CGAL/Mesh_3/Polyline_with_context.h \
    CGAL/Mesh_3/polylines_to_protect.h \
    CGAL/Mesh_3/Profile_counter.h \
    CGAL/Mesh_3/Profiling_tools.h \
    CGAL/Mesh_3/Protect_edges_sizing_field.h \
    CGAL/Mesh_3/radius_ratio.h \
    CGAL/Mesh_3/Refine_cells_3.h \
    CGAL/Mesh_3/Refine_facets_3.h \
    CGAL/Mesh_3/Refine_facets_manifold_base.h \
    CGAL/Mesh_3/Refine_tets_visitor.h \
    CGAL/Mesh_3/Robust_intersection_kernel.h \
    CGAL/Mesh_3/Robust_intersection_traits_3.h \
    CGAL/Mesh_3/search_for_connected_components_in_labeled_image.h \
    CGAL/Mesh_3/Sizing_grid.h \
    CGAL/Mesh_3/sliver_criteria.h \
    CGAL/Mesh_3/Sliver_perturber.h \
    CGAL/Mesh_3/Slivers_exuder.h \
    CGAL/Mesh_3/Slivers_exuder_cell_attributes_traits.h \
    CGAL/Mesh_3/squared_distance_Point_3_Triangle_3.h \
    CGAL/Mesh_3/tet_soup_to_c3t3.h \
    CGAL/Mesh_3/Triangle_accessor_primitive.h \
    CGAL/Mesh_3/Triangulation_helpers.h \
    CGAL/Mesh_3/Triangulation_sizing_field.h \
    CGAL/Mesh_3/Uniform_sizing_field.h \
    CGAL/Mesh_3/utilities.h \
    CGAL/Mesh_3/vertex_perturbation.h \
    CGAL/Mesh_3/Worksharing_data_structures.h \
    CGAL/Meshes/Double_map_container.h \
    CGAL/Meshes/Filtered_deque_container.h \
    CGAL/Meshes/Filtered_multimap_container.h \
    CGAL/Meshes/Filtered_queue_container.h \
    CGAL/Meshes/Simple_map_container.h \
    CGAL/Meshes/Simple_queue_container.h \
    CGAL/Meshes/Simple_set_container.h \
    CGAL/Meshes/Triangulation_mesher_level_traits_2.h \
    CGAL/Meshes/Triangulation_mesher_level_traits_3.h \
#    CGAL/Min_circle_2/Min_circle_2_adapterC2.h \
#    CGAL/Min_circle_2/Min_circle_2_adapterH2.h \
#    CGAL/Min_circle_2/Min_circle_2_impl.h \
#    CGAL/Min_circle_2/Optimisation_circle_2.h \
#    CGAL/Min_circle_2/Optimisation_circle_2_impl.h \
#    CGAL/Min_ellipse_2/Min_ellipse_2_adapterC2.h \
#    CGAL/Min_ellipse_2/Min_ellipse_2_adapterH2.h \
#    CGAL/Min_ellipse_2/Min_ellipse_2_impl.h \
#    CGAL/Min_ellipse_2/Optimisation_ellipse_2.h \
#    CGAL/Min_ellipse_2/Optimisation_ellipse_2_impl.h \
#    CGAL/Min_sphere_d/Min_sphere_d_impl.h \
#    CGAL/Min_sphere_d/Optimisation_sphere_d.h \
#    CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_configure.h \
#    CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_impl.h \
#    CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_pair.h \
#    CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_pivot_impl.h \
#    CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_support_set.h \
#    CGAL/Min_sphere_of_spheres_d/Min_sphere_of_spheres_d_support_set_impl.h \
#    CGAL/Minkowski_sum_2/AABB_collision_detector_2.h \
#    CGAL/Minkowski_sum_2/AABB_node_with_join.h \
#    CGAL/Minkowski_sum_2/AABB_segment_2_primitive.h \
#    CGAL/Minkowski_sum_2/AABB_traits_2.h \
#    CGAL/Minkowski_sum_2/AABB_traversal_traits_with_join.h \
#    CGAL/Minkowski_sum_2/AABB_tree_with_join.h \
#    CGAL/Minkowski_sum_2/Approx_offset_base_2.h \
#    CGAL/Minkowski_sum_2/Arr_labeled_traits_2.h \
#    CGAL/Minkowski_sum_2/Decomposition_strategy_adapter.h \
#    CGAL/Minkowski_sum_2/Exact_offset_base_2.h \
#    CGAL/Minkowski_sum_2/Hole_filter_2.h \
#    CGAL/Minkowski_sum_2/Labels.h \
#    CGAL/Minkowski_sum_2/Minkowski_sum_by_reduced_convolution_2.h \
#    CGAL/Minkowski_sum_2/Minkowski_sum_conv_2.h \
#    CGAL/Minkowski_sum_2/Minkowski_sum_decomp_2.h \
#    CGAL/Minkowski_sum_2/Offset_conv_2.h \
#    CGAL/Minkowski_sum_2/Offset_decomp_2.h \
#    CGAL/Minkowski_sum_2/Polygon_convex_decomposition.h \
#    CGAL/Minkowski_sum_2/Union_of_curve_cycles_2.h \
#    CGAL/Minkowski_sum_2/Union_of_cycles_2.h \
#    CGAL/Minkowski_sum_2/Union_of_segment_cycles_2.h \
#    CGAL/Minkowski_sum_3/bipartite_nary_union_sorted_combined.h \
#    CGAL/Minkowski_sum_3/Gaussian_map.h \
#    CGAL/Minkowski_sum_3/Gaussian_map_to_nef_3.h \
#    CGAL/Minkowski_sum_3/PointMark.h \
#    CGAL/Modular_arithmetic/Residue_type.h \
#    CGAL/Nef_2/Bounding_box_2.h \
#    CGAL/Nef_2/Constrained_triang_traits.h \
#    CGAL/Nef_2/debug.h \
#    CGAL/Nef_2/gen_point_location.h \
#    CGAL/Nef_2/geninfo.h \
#    CGAL/Nef_2/HDS_items.h \
#    CGAL/Nef_2/iterator_tools.h \
#    CGAL/Nef_2/Line_to_epoint.h \
#    CGAL/Nef_2/Object_handle.h \
#    CGAL/Nef_2/Object_index.h \
#    CGAL/Nef_2/PM_checker.h \
#    CGAL/Nef_2/PM_const_decorator.h \
#    CGAL/Nef_2/PM_decorator.h \
#    CGAL/Nef_2/PM_explorer.h \
#    CGAL/Nef_2/PM_io_parser.h \
#    CGAL/Nef_2/PM_overlayer.h \
#    CGAL/Nef_2/PM_persistent_PL.h \
#    CGAL/Nef_2/PM_point_locator.h \
#    CGAL/Nef_2/Polynomial.h \
#    CGAL/Nef_2/Polynomial_impl.h \
#    CGAL/Nef_2/Segment_overlay_traits.h \
#    CGAL/Nef_3/Binary_operation.h \
#    CGAL/Nef_3/binop_intersection_tests.h \
#    CGAL/Nef_3/bounded_side_3.h \
#    CGAL/Nef_3/Bounding_box_3.h \
#    CGAL/Nef_3/Combine_with_halfspace.h \
#    CGAL/Nef_3/Default_items.h \
#    CGAL/Nef_3/Edge_edge_overlay.h \
#    CGAL/Nef_3/Exact_triangulation_euclidean_traits_xy_3.h \
#    CGAL/Nef_3/Exact_triangulation_euclidean_traits_xz_3.h \
#    CGAL/Nef_3/Exact_triangulation_euclidean_traits_yz_3.h \
#    CGAL/Nef_3/Halfedge.h \
#    CGAL/Nef_3/Halffacet.h \
#    CGAL/Nef_3/ID_support_handler.h \
#    CGAL/Nef_3/Infimaximal_box.h \
#    CGAL/Nef_3/K3_tree.h \
#    CGAL/Nef_3/Mark_bounded_volumes.h \
#    CGAL/Nef_3/Nef_box.h \
#    CGAL/Nef_3/OGL_helper.h \
#    CGAL/Nef_3/Pluecker_line_3.h \
#    CGAL/Nef_3/polygon_mesh_to_nef_3.h \
#    CGAL/Nef_3/quotient_coordinates_to_homogeneous_point.h \
#    CGAL/Nef_3/SFace.h \
#    CGAL/Nef_3/SHalfedge.h \
#    CGAL/Nef_3/SHalfloop.h \
#    CGAL/Nef_3/shell_to_nef_3.h \
#    CGAL/Nef_3/SM_visualizor.h \
#    CGAL/Nef_3/SNC_const_decorator.h \
#    CGAL/Nef_3/SNC_constructor.h \
#    CGAL/Nef_3/SNC_decorator.h \
#    CGAL/Nef_3/SNC_decorator_traits.h \
#    CGAL/Nef_3/SNC_external_structure.h \
#    CGAL/Nef_3/SNC_FM_decorator.h \
#    CGAL/Nef_3/SNC_indexed_items.h \
#    CGAL/Nef_3/SNC_intersection.h \
#    CGAL/Nef_3/SNC_io_parser.h \
#    CGAL/Nef_3/SNC_items.h \
#    CGAL/Nef_3/SNC_iteration.h \
#    CGAL/Nef_3/SNC_k3_tree_traits.h \
#    CGAL/Nef_3/SNC_list.h \
#    CGAL/Nef_3/SNC_point_locator.h \
#    CGAL/Nef_3/SNC_simplify.h \
#    CGAL/Nef_3/SNC_SM_explorer.h \
#    CGAL/Nef_3/SNC_SM_overlayer.h \
#    CGAL/Nef_3/SNC_SM_visualizor.h \
#    CGAL/Nef_3/SNC_sphere_map.h \
#    CGAL/Nef_3/SNC_structure.h \
#    CGAL/Nef_3/Vertex.h \
#    CGAL/Nef_3/vertex_cycle_to_nef_3.h \
#    CGAL/Nef_3/Volume.h \
#    CGAL/Nef_S2/Generic_handle_map.h \
#    CGAL/Nef_S2/ID_support_handler.h \
#    CGAL/Nef_S2/leda_sphere_map.h \
#    CGAL/Nef_S2/Normalizing.h \
#    CGAL/Nef_S2/OGL_base_object.h \
#    CGAL/Nef_S2/SM_checker.h \
#    CGAL/Nef_S2/SM_const_decorator.h \
#    CGAL/Nef_S2/SM_constrained_triang_traits.h \
#    CGAL/Nef_S2/SM_decorator.h \
#    CGAL/Nef_S2/SM_decorator_traits.h \
#    CGAL/Nef_S2/SM_io_parser.h \
#    CGAL/Nef_S2/SM_items.h \
#    CGAL/Nef_S2/SM_iteration.h \
#    CGAL/Nef_S2/SM_list.h \
#    CGAL/Nef_S2/SM_overlayer.h \
#    CGAL/Nef_S2/SM_point_locator.h \
#    CGAL/Nef_S2/SM_triangulator.h \
#    CGAL/Nef_S2/SM_visualizor.h \
#    CGAL/Nef_S2/Sphere_circle.h \
#    CGAL/Nef_S2/Sphere_direction.h \
#    CGAL/Nef_S2/Sphere_geometry.h \
#    CGAL/Nef_S2/Sphere_geometry_OGL.h \
#    CGAL/Nef_S2/Sphere_map.h \
#    CGAL/Nef_S2/Sphere_point.h \
#    CGAL/Nef_S2/sphere_predicates.h \
#    CGAL/Nef_S2/Sphere_segment.h \
#    CGAL/Nef_S2/Sphere_triangle.h \
#    CGAL/NewKernel_d/LA_eigen/constructors.h \
#    CGAL/NewKernel_d/LA_eigen/LA.h \
#    CGAL/NewKernel_d/Types/Aff_transformation.h \
#    CGAL/NewKernel_d/Types/Hyperplane.h \
#    CGAL/NewKernel_d/Types/Iso_box.h \
#    CGAL/NewKernel_d/Types/Line.h \
#    CGAL/NewKernel_d/Types/Ray.h \
#    CGAL/NewKernel_d/Types/Segment.h \
#    CGAL/NewKernel_d/Types/Sphere.h \
#    CGAL/NewKernel_d/Types/Weighted_point.h \
#    CGAL/NewKernel_d/Vector/array.h \
#    CGAL/NewKernel_d/Vector/avx4.h \
#    CGAL/NewKernel_d/Vector/determinant_of_iterator_to_points_from_iterator_to_vectors.h \
#    CGAL/NewKernel_d/Vector/determinant_of_iterator_to_points_from_points.h \
#    CGAL/NewKernel_d/Vector/determinant_of_iterator_to_vectors_from_vectors.h \
#    CGAL/NewKernel_d/Vector/determinant_of_points_from_vectors.h \
#    CGAL/NewKernel_d/Vector/determinant_of_vectors_small_dim.h \
#    CGAL/NewKernel_d/Vector/determinant_of_vectors_small_dim_internal.h \
#    CGAL/NewKernel_d/Vector/mix.h \
#    CGAL/NewKernel_d/Vector/sse2.h \
#    CGAL/NewKernel_d/Vector/v2int.h \
#    CGAL/NewKernel_d/Vector/vector.h \
#    CGAL/NewKernel_d/Wrapper/Cartesian_wrap.h \
#    CGAL/NewKernel_d/Wrapper/Hyperplane_d.h \
#    CGAL/NewKernel_d/Wrapper/Point_d.h \
#    CGAL/NewKernel_d/Wrapper/Ref_count_obj.h \
#    CGAL/NewKernel_d/Wrapper/Segment_d.h \
#    CGAL/NewKernel_d/Wrapper/Sphere_d.h \
#    CGAL/NewKernel_d/Wrapper/Vector_d.h \
#    CGAL/NewKernel_d/Wrapper/Weighted_point_d.h \
#    CGAL/NewKernel_d/Cartesian_base.h \
#    CGAL/NewKernel_d/Cartesian_change_FT.h \
#    CGAL/NewKernel_d/Cartesian_complete.h \
#    CGAL/NewKernel_d/Cartesian_filter_K.h \
#    CGAL/NewKernel_d/Cartesian_filter_NT.h \
#    CGAL/NewKernel_d/Cartesian_LA_base.h \
#    CGAL/NewKernel_d/Cartesian_LA_functors.h \
#    CGAL/NewKernel_d/Cartesian_per_dimension.h \
#    CGAL/NewKernel_d/Cartesian_static_filters.h \
#    CGAL/NewKernel_d/Coaffine.h \
#    CGAL/NewKernel_d/Define_kernel_types.h \
#    CGAL/NewKernel_d/Dimension_base.h \
#    CGAL/NewKernel_d/Filtered_predicate2.h \
#    CGAL/NewKernel_d/function_objects_cartesian.h \
#    CGAL/NewKernel_d/functor_properties.h \
#    CGAL/NewKernel_d/functor_tags.h \
#    CGAL/NewKernel_d/Kernel_2_interface.h \
#    CGAL/NewKernel_d/Kernel_3_interface.h \
#    CGAL/NewKernel_d/Kernel_d_interface.h \
#    CGAL/NewKernel_d/Kernel_object_converter.h \
#    CGAL/NewKernel_d/KernelD_converter.h \
#    CGAL/NewKernel_d/Lazy_cartesian.h \
#    CGAL/NewKernel_d/static_int.h \
#    CGAL/NewKernel_d/store_kernel.h \
#    CGAL/NewKernel_d/utils.h \
#    CGAL/Number_types/internal_functions_comparison_root_of_2.h \
#    CGAL/OpenNL/bicgstab.h \
#    CGAL/OpenNL/blas.h \
#    CGAL/OpenNL/conjugate_gradient.h \
#    CGAL/OpenNL/full_vector.h \
#    CGAL/OpenNL/linear_solver.h \
#    CGAL/OpenNL/preconditioner.h \
#    CGAL/OpenNL/sparse_matrix.h \
#    CGAL/Optimisation/Access_coordinates_begin_2.h \
#    CGAL/Optimisation/Access_coordinates_begin_3.h \
#    CGAL/Optimisation/Access_coordinates_begin_d.h \
#    CGAL/Optimisation/Access_dimension_2.h \
#    CGAL/Optimisation/Access_dimension_3.h \
#    CGAL/Optimisation/Access_dimension_d.h \
#    CGAL/Optimisation/assertions.h \
#    CGAL/Optimisation/basic.h \
#    CGAL/Optimisation/Construct_point_2.h \
#    CGAL/Optimisation/Construct_point_3.h \
#    CGAL/Optimisation/Construct_point_d.h \
#    CGAL/Optimisation/debug.h \
#    CGAL/OTR_2/Cost.h \
#    CGAL/OTR_2/Reconstruction_edge_2.h \
#    CGAL/OTR_2/Reconstruction_face_base_2.h \
#    CGAL/OTR_2/Reconstruction_triangulation_2.h \
#    CGAL/OTR_2/Reconstruction_vertex_base_2.h \
#    CGAL/OTR_2/Sample.h \
#    CGAL/Partition_2/Circulator_pair.h \
#    CGAL/Partition_2/Indirect_edge_compare.h \
#    CGAL/Partition_2/Indirect_less_xy_2.h \
#    CGAL/Partition_2/Indirect_not_less_yx_2.h \
#    CGAL/Partition_2/is_degenerate_polygon_2.h \
#    CGAL/Partition_2/Iterator_list.h \
#    CGAL/Partition_2/Matrix.h \
#    CGAL/Partition_2/partition_approx_convex_2.h \
#    CGAL/Partition_2/partition_assertions.h \
#    CGAL/Partition_2/partition_greene_approx_convex_2.h \
#    CGAL/Partition_2/Partition_opt_cvx_diagonal_list.h \
#    CGAL/Partition_2/Partition_opt_cvx_edge.h \
#    CGAL/Partition_2/Partition_opt_cvx_vertex.h \
#    CGAL/Partition_2/partition_optimal_convex_2.h \
#    CGAL/Partition_2/Partition_traits_2_base.h \
#    CGAL/Partition_2/Partition_vertex_map.h \
#    CGAL/Partition_2/partition_y_monotone_2.h \
#    CGAL/Partition_2/Partitioned_polygon_2.h \
#    CGAL/Partition_2/Point_pair_less_xy_2.h \
#    CGAL/Partition_2/Rotation_tree_2.h \
#    CGAL/Partition_2/Rotation_tree_2_impl.h \
#    CGAL/Partition_2/Rotation_tree_node_2.h \
#    CGAL/Partition_2/Segment_less_yx_2.h \
#    CGAL/Partition_2/Triangulation_indirect_traits_2.h \
#    CGAL/Partition_2/Turn_reverser.h \
#    CGAL/Partition_2/Vertex_visibility_graph_2.h \
#    CGAL/Partition_2/Vertex_visibility_graph_2_impl.h \
#    CGAL/Periodic_3_mesh_3/IO/File_medit.h \
#    CGAL/Periodic_3_mesh_3/config.h \
#    CGAL/Periodic_3_mesh_3/Protect_edges_sizing_field.h \
#    CGAL/Point_set_3/IO.h \
#    CGAL/Point_set_processing_3/internal/Voronoi_covariance_3/voronoi_covariance_3.h \
#    CGAL/Point_set_processing_3/internal/Voronoi_covariance_3/voronoi_covariance_sphere_3.h \
#    CGAL/Point_set_processing_3/internal/Rich_grid.h \
#    CGAL/Polygon_2/Polygon_2_algorithms_impl.h \
#    CGAL/Polygon_2/Polygon_2_edge_circulator.h \
#    CGAL/Polygon_2/Polygon_2_edge_iterator.h \
#    CGAL/Polygon_2/Polygon_2_impl.h \
#    CGAL/Polygon_2/Polygon_2_simplicity.h \
#    CGAL/Polygon_2/Polygon_2_vertex_circulator.h \
#    CGAL/Polygon_2/polygon_assertions.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/Face_graph_output_builder.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/face_graph_utils.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/intersect_triangle_and_segment_3.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_callbacks.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_impl.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_nodes.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_of_coplanar_triangles_3.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/Intersection_type.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/Output_builder_for_autorefinement.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/predicates.h \
    CGAL/Polygon_mesh_processing/internal/Corefinement/Visitor.h \
    CGAL/Polygon_mesh_processing/internal/Hole_filling/experimental/experimental_code.h \
    CGAL/Polygon_mesh_processing/internal/Hole_filling/do_not_use_DT3.h \
    CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polygon_mesh.h \
    CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polyline.h \
    CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/AABB_filtered_projection_traits.h \
    CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/remesh_impl.h \
    CGAL/Polygon_mesh_processing/internal/Polygon_mesh_slicer/Axis_parallel_plane_traits.h \
    CGAL/Polygon_mesh_processing/internal/Polygon_mesh_slicer/Traversal_traits.h \
    CGAL/Polygon_mesh_processing/internal/Side_of_triangle_mesh/Point_inside_vertical_ray_cast.h \
    CGAL/Polygon_mesh_processing/internal/Side_of_triangle_mesh/Ray_3_Triangle_3_traversal_traits.h \
    CGAL/Polygon_mesh_processing/internal/do_no_use_CDT2.h \
    CGAL/Polygon_mesh_processing/internal/fair_impl.h \
    CGAL/Polygon_mesh_processing/internal/mesh_to_point_set_hausdorff_distance.h \
    CGAL/Polygon_mesh_processing/internal/named_function_params.h \
    CGAL/Polygon_mesh_processing/internal/named_params_helper.h \
    CGAL/Polygon_mesh_processing/internal/refine_impl.h \
    CGAL/Polygon_mesh_processing/internal/repair_extra.h \
    CGAL/Polygon_mesh_processing/bbox.h \
    CGAL/Polygon_mesh_processing/border.h \
    CGAL/Polygon_mesh_processing/clip.h \
    CGAL/Polygon_mesh_processing/compute_normal.h \
    CGAL/Polygon_mesh_processing/connected_components.h \
    CGAL/Polygon_mesh_processing/corefinement.h \
    CGAL/Polygon_mesh_processing/detect_features.h \
    CGAL/Polygon_mesh_processing/distance.h \
    CGAL/Polygon_mesh_processing/extrude.h \
    CGAL/Polygon_mesh_processing/fair.h \
    CGAL/Polygon_mesh_processing/intersection.h \
    CGAL/Polygon_mesh_processing/measure.h \
    CGAL/Polygon_mesh_processing/orient_polygon_soup.h \
    CGAL/Polygon_mesh_processing/orientation.h \
    CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h \
    CGAL/Polygon_mesh_processing/random_perturbation.h \
    CGAL/Polygon_mesh_processing/refine.h \
    CGAL/Polygon_mesh_processing/remesh.h \
    CGAL/Polygon_mesh_processing/repair.h \
    CGAL/Polygon_mesh_processing/self_intersections.h \
    CGAL/Polygon_mesh_processing/stitch_borders.h \
    CGAL/Polygon_mesh_processing/transform.h \
    CGAL/Polygon_mesh_processing/triangulate_faces.h \
    CGAL/Polygon_mesh_processing/triangulate_hole.h \
    CGAL/Polygon_mesh_processing/Weights.h \
    CGAL/Polyline_simplification_2/Hybrid_squared_distance_cost.h \
    CGAL/Polyline_simplification_2/Scaled_squared_distance_cost.h \
    CGAL/Polyline_simplification_2/simplify.h \
    CGAL/Polyline_simplification_2/Squared_distance_cost.h \
    CGAL/Polyline_simplification_2/Stop_above_cost_threshold.h \
    CGAL/Polyline_simplification_2/Stop_below_count_ratio_threshold.h \
    CGAL/Polyline_simplification_2/Stop_below_count_threshold.h \
    CGAL/Polyline_simplification_2/Vertex_base_2.h \
#    CGAL/Polynomial/Algebraic_structure_traits.h \
#    CGAL/Polynomial/bezout_matrix.h \
#    CGAL/Polynomial/Cached_extended_euclidean_algorithm.h \
#    CGAL/Polynomial/Chinese_remainder_traits.h \
#    CGAL/Polynomial/Coercion_traits.h \
#    CGAL/Polynomial/Degree.h \
#    CGAL/Polynomial/determinant.h \
#    CGAL/Polynomial/Fraction_traits.h \
#    CGAL/Polynomial/fwd.h \
#    CGAL/Polynomial/Get_arithmetic_kernel.h \
#    CGAL/Polynomial/hgdelta_update.h \
#    CGAL/Polynomial/Interpolator.h \
#    CGAL/Polynomial/misc.h \
#    CGAL/Polynomial/modular_filter.h \
#    CGAL/Polynomial/modular_gcd.h \
#    CGAL/Polynomial/modular_gcd_utcf_algorithm_M.h \
#    CGAL/Polynomial/modular_gcd_utcf_dfai.h \
#    CGAL/Polynomial/modular_gcd_utils.h \
#    CGAL/Polynomial/Modular_traits.h \
#    CGAL/Polynomial/Monomial_representation.h \
#    CGAL/Polynomial/polynomial_gcd.h \
#    CGAL/Polynomial/polynomial_gcd_implementations.h \
#    CGAL/Polynomial/polynomial_gcd_ntl.h \
#    CGAL/Polynomial/Polynomial_type.h \
#    CGAL/Polynomial/prs_resultant.h \
#    CGAL/Polynomial/Real_embeddable_traits.h \
#    CGAL/Polynomial/resultant.h \
#    CGAL/Polynomial/Scalar_factor_traits.h \
#    CGAL/Polynomial/square_free_factorize.h \
#    CGAL/Polynomial/sturm_habicht_sequence.h \
#    CGAL/Polynomial/subresultants.h \
#    CGAL/predicates/kernel_ftC2.h \
#    CGAL/predicates/kernel_ftC3.h \
#    CGAL/predicates/Polygon_offset_pred_ftC2.h \
#    CGAL/predicates/predicates_for_mixed_complex_3.h \
#    CGAL/predicates/predicates_for_voronoi_intersection_cartesian_2_3.h \
#    CGAL/predicates/sign_of_determinant.h \
#    CGAL/predicates/Straight_skeleton_pred_ftC2.h \
#    CGAL/QP_solver/assertions.h \
#    CGAL/QP_solver/basic.h \
#    CGAL/QP_solver/debug.h \
#    CGAL/QP_solver/functors.h \
#    CGAL/QP_solver/Initialization.h \
#    CGAL/QP_solver/QP__filtered_base.h \
#    CGAL/QP_solver/QP__filtered_base_impl.h \
#    CGAL/QP_solver/QP__partial_base.h \
#    CGAL/QP_solver/QP_basis_inverse.h \
#    CGAL/QP_solver/QP_basis_inverse_impl.h \
#    CGAL/QP_solver/QP_exact_bland_pricing.h \
#    CGAL/QP_solver/QP_full_exact_pricing.h \
#    CGAL/QP_solver/QP_full_filtered_pricing.h \
#    CGAL/QP_solver/QP_functions_impl.h \
#    CGAL/QP_solver/QP_partial_exact_pricing.h \
#    CGAL/QP_solver/QP_partial_filtered_pricing.h \
#    CGAL/QP_solver/QP_pricing_strategy.h \
#    CGAL/QP_solver/QP_solution_impl.h \
#    CGAL/QP_solver/QP_solver.h \
#    CGAL/QP_solver/QP_solver_bounds_impl.h \
#    CGAL/QP_solver/QP_solver_impl.h \
#    CGAL/QP_solver/QP_solver_nonstandardform_impl.h \
#    CGAL/QP_solver/Unbounded_direction.h \
#    CGAL/Qt/AlphaShapeGraphicsItem.h \
#    CGAL/Qt/ApolloniusGraphGraphicsItem.h \
#    CGAL/Qt/Basic_viewer_qt.h \
#    CGAL/Qt/camera.h \
#    CGAL/Qt/camera_impl.h \
#    CGAL/Qt/CGAL_Qt_config.h \
#    CGAL/Qt/CircularArcGraphicsItem.h \
#    CGAL/Qt/ConstrainedTriangulationGraphicsItem.h \
#    CGAL/Qt/constraint.h \
#    CGAL/Qt/constraint_impl.h \
#    CGAL/Qt/Converter.h \
#    CGAL/Qt/CreateOpenGLContext.h \
#    CGAL/Qt/debug.h \
#    CGAL/Qt/debug_impl.h \
#    CGAL/Qt/DelaunayMeshTriangulationGraphicsItem.h \
#    CGAL/Qt/DemosMainWindow.h \
#    CGAL/Qt/DemosMainWindow_impl.h \
#    CGAL/Qt/domUtils.h \
#    CGAL/Qt/frame.h \
#    CGAL/Qt/frame_impl.h \
#    CGAL/Qt/GraphicsItem.h \
#    CGAL/Qt/GraphicsViewCircleInput.h \
#    CGAL/Qt/GraphicsViewCircularArcInput.h \
#    CGAL/Qt/GraphicsViewInput.h \
#    CGAL/Qt/GraphicsViewIsoRectangleInput.h \
#    CGAL/Qt/GraphicsViewLineInput.h \
#    CGAL/Qt/GraphicsViewNavigation.h \
#    CGAL/Qt/GraphicsViewNavigation_impl.h \
#    CGAL/Qt/GraphicsViewPointInput.h \
#    CGAL/Qt/GraphicsViewPolygonWithHolesInput.h \
#    CGAL/Qt/GraphicsViewPolylineInput.h \
#    CGAL/Qt/GraphicsViewPolylineInput_impl.h \
#    CGAL/Qt/image_interface.h \
#    CGAL/Qt/keyFrameInterpolator.h \
#    CGAL/Qt/keyFrameInterpolator_impl.h \
#    CGAL/Qt/LineGraphicsItem.h \
#    CGAL/Qt/manipulatedCameraFrame.h \
#    CGAL/Qt/manipulatedCameraFrame_impl.h \
#    CGAL/Qt/manipulatedFrame.h \
#    CGAL/Qt/manipulatedFrame_impl.h \
#    CGAL/Qt/mouseGrabber.h \
#    CGAL/Qt/mouseGrabber_impl.h \
#    CGAL/Qt/PainterOstream.h \
#    CGAL/Qt/PointsGraphicsItem.h \
#    CGAL/Qt/PointsInKdTreeGraphicsItem.h \
#    CGAL/Qt/PolygonGraphicsItem.h \
#    CGAL/Qt/PolygonWithHolesGraphicsItem.h \
#    CGAL/Qt/PolylinesGraphicsItem.h \
#    CGAL/Qt/PowerdiagramGraphicsItem.h \
#    CGAL/Qt/qglviewer.h \
#    CGAL/Qt/qglviewer_impl.h \
#    CGAL/Qt/qglviewer_impl_list.h \
#    CGAL/Qt/quaternion.h \
#    CGAL/Qt/quaternion_impl.h \
#    CGAL/Qt/RegularGridGraphicsItem.h \
#    CGAL/Qt/RegularGridVectorFieldGraphicsItem.h \
#    CGAL/Qt/RegularTriangulationGraphicsItem.h \
#    CGAL/Qt/resources.h \
#    CGAL/Qt/resources_impl.h \
#    CGAL/Qt/SegmentDelaunayGraphGraphicsItem.h \
#    CGAL/Qt/SegmentDelaunayGraphLinfGraphicsItem.h \
#    CGAL/Qt/SegmentsGraphicsItem.h \
#    CGAL/Qt/StreamLinesGraphicsItem.h \
#    CGAL/Qt/TriangulationGraphicsItem.h \
#    CGAL/Qt/utility.h \
#    CGAL/Qt/utility_impl.h \
#    CGAL/Qt/vec.h \
#    CGAL/Qt/vec_impl.h \
#    CGAL/Qt/viewer_actions.h \
#    CGAL/Qt/VoronoiGraphicsItem.h \
#    CGAL/RS/ak_1.h \
#    CGAL/RS/ak_z_1.h \
#    CGAL/RS/algebraic_1.h \
#    CGAL/RS/algebraic_z_1.h \
#    CGAL/RS/bisection_refiner_1.h \
#    CGAL/RS/comparator_1.h \
#    CGAL/RS/dyadic.h \
#    CGAL/RS/exact_signat_1.h \
#    CGAL/RS/functors_1.h \
#    CGAL/RS/functors_z_1.h \
#    CGAL/RS/Gmpfr_make_unique.h \
#    CGAL/RS/polynomial_converter_1.h \
#    CGAL/RS/rs23_k_isolator_1.h \
#    CGAL/RS/rs2_calls.h \
#    CGAL/RS/rs2_isolator_1.h \
#    CGAL/RS/rs3_k_refiner_1.h \
#    CGAL/RS/rs3_refiner_1.h \
#    CGAL/RS/signat_1.h \
#    CGAL/Scale_space_reconstruction_3/internal/Auto_count.h \
#    CGAL/Scale_space_reconstruction_3/Advancing_front_mesher.h \
#    CGAL/Scale_space_reconstruction_3/Alpha_shape_mesher.h \
#    CGAL/Scale_space_reconstruction_3/Jet_smoother.h \
#    CGAL/Scale_space_reconstruction_3/Shape_construction_3.h \
#    CGAL/Scale_space_reconstruction_3/Weighted_PCA_smoother.h \
#    CGAL/Segment_Delaunay_graph_2/Are_parallel_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Are_same_points_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Are_same_segments_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Arrangement_enum.h \
#    CGAL/Segment_Delaunay_graph_2/Arrangement_type_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Arrangement_type_non_intersecting_C2.h \
#    CGAL/Segment_Delaunay_graph_2/basic.h \
#    CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Cartesian_converter.h \
#    CGAL/Segment_Delaunay_graph_2/Compare_x_2.h \
#    CGAL/Segment_Delaunay_graph_2/Compare_y_2.h \
#    CGAL/Segment_Delaunay_graph_2/Construct_storage_site_2.h \
#    CGAL/Segment_Delaunay_graph_2/Construct_storage_site_with_info_2.h \
#    CGAL/Segment_Delaunay_graph_2/Constructions_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Filtered_traits_base_2.h \
#    CGAL/Segment_Delaunay_graph_2/Filtered_traits_concept_check_tags.h \
#    CGAL/Segment_Delaunay_graph_2/Finite_edge_interior_conflict_C2.h \
#    CGAL/Segment_Delaunay_graph_2/in_place_edge_list.h \
#    CGAL/Segment_Delaunay_graph_2/Infinite_edge_interior_conflict_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Is_degenerate_edge_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Kernel_wrapper_2.h \
#    CGAL/Segment_Delaunay_graph_2/Orientation_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Oriented_side_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Oriented_side_of_bisector_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Predicates_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Segment_Delaunay_graph_2_impl.h \
#    CGAL/Segment_Delaunay_graph_2/Segment_Delaunay_graph_hierarchy_2_impl.h \
#    CGAL/Segment_Delaunay_graph_2/Sqrt_extension_2.h \
#    CGAL/Segment_Delaunay_graph_2/Traits_base_2.h \
#    CGAL/Segment_Delaunay_graph_2/Traits_wrapper_2.h \
#    CGAL/Segment_Delaunay_graph_2/Triangulation_face_base_with_edges_2.h \
#    CGAL/Segment_Delaunay_graph_2/Vertex_conflict_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_ring_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_sqrt_field_C2.h \
#    CGAL/Segment_Delaunay_graph_2/Voronoi_vertex_sqrt_field_new_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/basic.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Basic_predicates_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Bisector_Linf.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Constructions_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Filtered_traits_base_2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Finite_edge_interior_conflict_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Infinite_edge_interior_conflict_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Orientation_Linf_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Oriented_side_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Oriented_side_of_bisector_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Predicates_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Segment_Delaunay_graph_Linf_2_impl.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Segment_Delaunay_graph_Linf_hierarchy_2_impl.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Traits_base_2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Vertex_conflict_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_ring_C2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2/Voronoi_vertex_sqrt_field_new_C2.h \
#    CGAL/Set_movable_separability_2/internal/Circle_arrangment.h \
#    CGAL/Set_movable_separability_2/internal/Utils.h \
#    CGAL/Set_movable_separability_2/Single_mold_translational_casting/is_pullout_direction.h \
#    CGAL/Set_movable_separability_2/Single_mold_translational_casting/pullout_directions.h \
#    CGAL/Set_movable_separability_2/Single_mold_translational_casting/top_edges.h \
#    CGAL/Shape_detection_3/Cone.h \
#    CGAL/Shape_detection_3/Cylinder.h \
#    CGAL/Shape_detection_3/Efficient_RANSAC.h \
#    CGAL/Shape_detection_3/Efficient_RANSAC_traits.h \
#    CGAL/Shape_detection_3/Octree.h \
#    CGAL/Shape_detection_3/Plane.h \
#    CGAL/Shape_detection_3/property_maps.h \
#    CGAL/Shape_detection_3/Region_growing.h \
#    CGAL/Shape_detection_3/Shape_base.h \
#    CGAL/Shape_detection_3/Shape_detection_traits.h \
#    CGAL/Shape_detection_3/Sphere.h \
#    CGAL/Shape_detection_3/Torus.h \
#    CGAL/Sqrt_extension/Algebraic_extension_traits.h \
#    CGAL/Sqrt_extension/Algebraic_structure_traits.h \
#    CGAL/Sqrt_extension/Chinese_remainder_traits.h \
#    CGAL/Sqrt_extension/Coercion_traits.h \
#    CGAL/Sqrt_extension/convert_to_bfi.h \
#    CGAL/Sqrt_extension/Eigen_NumTraits.h \
#    CGAL/Sqrt_extension/Fraction_traits.h \
#    CGAL/Sqrt_extension/Get_arithmetic_kernel.h \
#    CGAL/Sqrt_extension/io.h \
#    CGAL/Sqrt_extension/Modular_traits.h \
#    CGAL/Sqrt_extension/Real_embeddable_traits.h \
#    CGAL/Sqrt_extension/Scalar_factor_traits.h \
#    CGAL/Sqrt_extension/Sqrt_extension_type.h \
#    CGAL/Sqrt_extension/Wang_traits.h \
#    CGAL/Straight_skeleton_2/assertions.h \
#    CGAL/Straight_skeleton_2/debug.h \
#    CGAL/Straight_skeleton_2/Polygon_offset_builder_2_impl.h \
#    CGAL/Straight_skeleton_2/Straight_skeleton_aux.h \
#    CGAL/Straight_skeleton_2/Straight_skeleton_builder_2_impl.h \
#    CGAL/Straight_skeleton_2/Straight_skeleton_builder_events_2.h \
#    CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h \
#    CGAL/Straight_skeleton_2/test.h \
#    CGAL/Subdivision_method_3/internal/Euler_extensions.h \
#    CGAL/Subdivision_method_3/internal/subdivision_hosts_impl_3.h \
#    CGAL/Subdivision_method_3/subdivision_hosts_3.h \
#    CGAL/Subdivision_method_3/subdivision_masks_3.h \
#    CGAL/Subdivision_method_3/subdivision_methods_3.h \
    CGAL/Surface_mesh/IO.h \
    CGAL/Surface_mesh/Properties.h \
    CGAL/Surface_mesh/Surface_mesh.h \
    CGAL/Surface_mesh/Surface_mesh_fwd.h \
    CGAL/Surface_mesh_parameterization/internal/angles.h \
    CGAL/Surface_mesh_parameterization/internal/Bool_property_map.h \
    CGAL/Surface_mesh_parameterization/internal/Containers_filler.h \
    CGAL/Surface_mesh_parameterization/internal/kernel_traits.h \
    CGAL/Surface_mesh_parameterization/internal/orbifold_cone_helper.h \
    CGAL/Surface_mesh_parameterization/internal/validity.h \
    CGAL/Surface_mesh_parameterization/IO/File_off.h \
    CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h \
    CGAL/Surface_mesh_parameterization/Barycentric_mapping_parameterizer_3.h \
    CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h \
    CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h \
    CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h \
    CGAL/Surface_mesh_parameterization/Error_code.h \
    CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h \
    CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h \
    CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h \
    CGAL/Surface_mesh_parameterization/MVC_post_processor_3.h \
    CGAL/Surface_mesh_parameterization/orbifold_enums.h \
    CGAL/Surface_mesh_parameterization/orbifold_shortest_path.h \
    CGAL/Surface_mesh_parameterization/Orbifold_Tutte_parameterizer_3.h \
    CGAL/Surface_mesh_parameterization/parameterize.h \
    CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h \
    CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h \
    CGAL/Surface_mesh_shortest_path/internal/Cone_expansion_event.h \
    CGAL/Surface_mesh_shortest_path/internal/Cone_tree.h \
    CGAL/Surface_mesh_shortest_path/internal/misc_functions.h \
    CGAL/Surface_mesh_shortest_path/barycentric.h \
    CGAL/Surface_mesh_shortest_path/function_objects.h \
    CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path.h \
    CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path_traits.h \
    CGAL/Surface_mesh_simplification/Detail/Common.h \
    CGAL/Surface_mesh_simplification/Detail/Edge_collapse.h \
    CGAL/Surface_mesh_simplification/Detail/Edge_collapse_impl.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core_impl.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile_impl.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_params.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h \
    CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h \
    CGAL/Surface_mesh_simplification/edge_collapse.h \
    CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h \
    CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h \
#    CGAL/Surface_mesher/Combining_oracle.h \
#    CGAL/Surface_mesher/Has_edges.h \
#    CGAL/Surface_mesher/Implicit_surface_oracle_3.h \
#    CGAL/Surface_mesher/Intersection_data_structure_3.h \
#    CGAL/Surface_mesher/Null_oracle_visitor.h \
#    CGAL/Surface_mesher/Point_surface_indices_oracle_visitor.h \
#    CGAL/Surface_mesher/Poisson_implicit_surface_oracle_3.h \
#    CGAL/Surface_mesher/Polyhedral_oracle.h \
#    CGAL/Surface_mesher/Profile_counter.h \
#    CGAL/Surface_mesher/Profile_timer.h \
#    CGAL/Surface_mesher/Sphere_oracle_3.h \
#    CGAL/Surface_mesher/Standard_criteria.h \
#    CGAL/Surface_mesher/Surface_mesher.h \
#    CGAL/Surface_mesher/Surface_mesher_edges_level.h \
#    CGAL/Surface_mesher/Surface_mesher_edges_level_visitor.h \
#    CGAL/Surface_mesher/Surface_mesher_manifold.h \
#    CGAL/Surface_mesher/Surface_mesher_regular_edges.h \
#    CGAL/Surface_mesher/Surface_mesher_visitor.h \
#    CGAL/Surface_mesher/Types_generators.h \
#    CGAL/Surface_mesher/Verbose_flag.h \
#    CGAL/Surface_mesher/Vertices_on_the_same_psc_element_criterion.h \
#    CGAL/Surface_mesher/Vertices_on_the_same_surface_criterion.h \
#    CGAL/Surface_sweep_2/Arr_basic_insertion_traits_2.h \
#    CGAL/Surface_sweep_2/Arr_batched_pl_ss_visitor.h \
#    CGAL/Surface_sweep_2/Arr_construction_event.h \
#    CGAL/Surface_sweep_2/Arr_construction_event_base.h \
#    CGAL/Surface_sweep_2/Arr_construction_ss_visitor.h \
#    CGAL/Surface_sweep_2/Arr_construction_subcurve.h \
#    CGAL/Surface_sweep_2/Arr_default_overlay_traits_base.h \
#    CGAL/Surface_sweep_2/Arr_insertion_ss_visitor.h \
#    CGAL/Surface_sweep_2/Arr_insertion_traits_2.h \
#    CGAL/Surface_sweep_2/Arr_no_intersection_insertion_ss_visitor.h \
#    CGAL/Surface_sweep_2/Arr_overlay_event.h \
#    CGAL/Surface_sweep_2/Arr_overlay_ss_visitor.h \
#    CGAL/Surface_sweep_2/Arr_overlay_subcurve.h \
#    CGAL/Surface_sweep_2/Arr_overlay_traits_2.h \
#    CGAL/Surface_sweep_2/Arr_vert_decomp_ss_visitor.h \
#    CGAL/Surface_sweep_2/Curve_comparer.h \
#    CGAL/Surface_sweep_2/Curve_pair.h \
#    CGAL/Surface_sweep_2/Default_event.h \
#    CGAL/Surface_sweep_2/Default_event_base.h \
#    CGAL/Surface_sweep_2/Default_subcurve.h \
#    CGAL/Surface_sweep_2/Default_visitor.h \
#    CGAL/Surface_sweep_2/Default_visitor_base.h \
#    CGAL/Surface_sweep_2/Do_interior_intersect_visitor.h \
#    CGAL/Surface_sweep_2/Event_comparer.h \
#    CGAL/Surface_sweep_2/Intersection_points_visitor.h \
#    CGAL/Surface_sweep_2/No_intersection_surface_sweep_2_impl.h \
#    CGAL/Surface_sweep_2/No_overlap_event.h \
#    CGAL/Surface_sweep_2/No_overlap_event_base.h \
#    CGAL/Surface_sweep_2/No_overlap_subcurve.h \
#    CGAL/Surface_sweep_2/Subcurves_visitor.h \
#    CGAL/Surface_sweep_2/Surface_sweep_2_debug.h \
#    CGAL/Surface_sweep_2/Surface_sweep_2_impl.h \
#    CGAL/Surface_sweep_2/Surface_sweep_2_utils.h \
#    CGAL/Three/Buffer_objects.h \
#    CGAL/Three/Edge_container.h \
#    CGAL/Three/exceptions.h \
#    CGAL/Three/Point_container.h \
#    CGAL/Three/Polyhedron_demo_io_plugin_interface.h \
#    CGAL/Three/Polyhedron_demo_plugin_helper.h \
#    CGAL/Three/Polyhedron_demo_plugin_interface.h \
#    CGAL/Three/Primitive_container.h \
#    CGAL/Three/Scene_draw_interface.h \
#    CGAL/Three/Scene_group_item.h \
#    CGAL/Three/Scene_interface.h \
#    CGAL/Three/Scene_item.h \
#    CGAL/Three/Scene_item_config.h \
#    CGAL/Three/Scene_item_rendering_helper.h \
#    CGAL/Three/Scene_item_with_properties.h \
#    CGAL/Three/Scene_print_item_interface.h \
#    CGAL/Three/Scene_transparent_interface.h \
#    CGAL/Three/Scene_zoomable_item_interface.h \
#    CGAL/Three/TextRenderer.h \
#    CGAL/Three/Three.h \
#    CGAL/Three/Triangle_container.h \
#    CGAL/Three/Viewer_config.h \
#    CGAL/Three/Viewer_interface.h \
#    CGAL/Tools/chained_map.h \
#    CGAL/Triangulation_2/insert_constraints.h \
#    CGAL/Visibility_2/visibility_utils.h \
    CGAL/Voronoi_diagram_2/Accessor.h \
    CGAL/Voronoi_diagram_2/Adaptation_traits_base_2.h \
    CGAL/Voronoi_diagram_2/Adaptation_traits_functors.h \
    CGAL/Voronoi_diagram_2/Apollonius_graph_degeneracy_testers.h \
    CGAL/Voronoi_diagram_2/Apollonius_graph_nearest_site_2.h \
    CGAL/Voronoi_diagram_2/basic.h \
    CGAL/Voronoi_diagram_2/Cached_degeneracy_testers.h \
    CGAL/Voronoi_diagram_2/Circulator_adaptors.h \
    CGAL/Voronoi_diagram_2/Connected_components.h \
    CGAL/Voronoi_diagram_2/Construct_dual_points.h \
    CGAL/Voronoi_diagram_2/Default_site_inserters.h \
    CGAL/Voronoi_diagram_2/Default_site_removers.h \
    CGAL/Voronoi_diagram_2/Degeneracy_tester_binders.h \
    CGAL/Voronoi_diagram_2/Delaunay_triangulation_degeneracy_testers.h \
    CGAL/Voronoi_diagram_2/Delaunay_triangulation_nearest_site_2.h \
    CGAL/Voronoi_diagram_2/Dummy_iterator.h \
    CGAL/Voronoi_diagram_2/Edge_less.h \
    CGAL/Voronoi_diagram_2/Face.h \
    CGAL/Voronoi_diagram_2/Finder_classes.h \
    CGAL/Voronoi_diagram_2/Halfedge.h \
    CGAL/Voronoi_diagram_2/Handle_adaptor.h \
    CGAL/Voronoi_diagram_2/Identity_rejectors.h \
    CGAL/Voronoi_diagram_2/Iterator_adaptors.h \
    CGAL/Voronoi_diagram_2/Policy_base.h \
    CGAL/Voronoi_diagram_2/Regular_triangulation_degeneracy_testers.h \
    CGAL/Voronoi_diagram_2/Regular_triangulation_nearest_site_2.h \
    CGAL/Voronoi_diagram_2/Segment_Delaunay_graph_degeneracy_testers.h \
    CGAL/Voronoi_diagram_2/Segment_Delaunay_graph_nearest_site_2.h \
    CGAL/Voronoi_diagram_2/Site_accessors.h \
    CGAL/Voronoi_diagram_2/Unbounded_edges.h \
    CGAL/Voronoi_diagram_2/Unbounded_faces.h \
    CGAL/Voronoi_diagram_2/Validity_testers.h \
    CGAL/Voronoi_diagram_2/Vertex.h \
#    CGAL/AABB_face_graph_triangle_primitive.h \
#    CGAL/AABB_halfedge_graph_segment_primitive.h \
#    CGAL/AABB_polyhedral_oracle.h \
#    CGAL/AABB_polyhedron_segment_primitive.h \
#    CGAL/AABB_polyhedron_triangle_primitive.h \
#    CGAL/AABB_primitive.h \
#    CGAL/AABB_segment_primitive.h \
#    CGAL/AABB_traits.h \
#    CGAL/AABB_tree.h \
#    CGAL/AABB_triangle_primitive.h \
#    CGAL/AABB_triangulation_3_triangle_primitive.h \
#    CGAL/Advancing_front_surface_reconstruction.h \
#    CGAL/Advancing_front_surface_reconstruction_cell_base_3.h \
#    CGAL/Advancing_front_surface_reconstruction_vertex_base_3.h \
#    CGAL/Aff_transformation_2.h \
#    CGAL/Aff_transformation_3.h \
#    CGAL/aff_transformation_tags.h \
#    CGAL/aff_transformation_tags_impl.h \
#    CGAL/Algebraic_extension_traits.h \
#    CGAL/Algebraic_kernel_converter.h \
#    CGAL/Algebraic_kernel_d_1.h \
#    CGAL/Algebraic_kernel_d_2.h \
#    CGAL/Algebraic_kernel_for_circles_2_2.h \
#    CGAL/Algebraic_kernel_for_spheres_2_3.h \
#    CGAL/Algebraic_kernel_rs_gmpq_d_1.h \
#    CGAL/Algebraic_kernel_rs_gmpz_d_1.h \
#    CGAL/Algebraic_structure_traits.h \
#    CGAL/algorithm.h \
#    CGAL/all_furthest_neighbors_2.h \
#    CGAL/Alpha_shape_2.h \
#    CGAL/Alpha_shape_3.h \
#    CGAL/Alpha_shape_cell_base_3.h \
#    CGAL/Alpha_shape_euclidean_traits_2.h \
#    CGAL/Alpha_shape_euclidean_traits_3.h \
#    CGAL/Alpha_shape_face_base_2.h \
#    CGAL/Alpha_shape_vertex_base_2.h \
#    CGAL/Alpha_shape_vertex_base_3.h \
#    CGAL/Apollonius_graph_2.h \
#    CGAL/Apollonius_graph_adaptation_policies_2.h \
#    CGAL/Apollonius_graph_adaptation_traits_2.h \
#    CGAL/Apollonius_graph_data_structure_2.h \
#    CGAL/Apollonius_graph_filtered_traits_2.h \
#    CGAL/Apollonius_graph_hierarchy_2.h \
#    CGAL/Apollonius_graph_hierarchy_vertex_base_2.h \
#    CGAL/Apollonius_graph_traits_2.h \
#    CGAL/Apollonius_graph_vertex_base_2.h \
#    CGAL/Apollonius_site_2.h \
#    CGAL/apply_to_range.h \
#    CGAL/Approximate_min_ellipsoid_d.h \
#    CGAL/Approximate_min_ellipsoid_d_traits_2.h \
#    CGAL/Approximate_min_ellipsoid_d_traits_3.h \
#    CGAL/Approximate_min_ellipsoid_d_traits_d.h \
#    CGAL/approximated_offset_2.h \
#    CGAL/argument_swaps.h \
#    CGAL/Arithmetic_kernel.h \
#    CGAL/Arr_accessor.h \
#    CGAL/Arr_algebraic_segment_traits_2.h \
#    CGAL/Arr_batched_point_location.h \
#    CGAL/Arr_Bezier_curve_traits_2.h \
#    CGAL/Arr_bounded_planar_topology_traits_2.h \
#    CGAL/Arr_circle_segment_traits_2.h \
#    CGAL/Arr_circular_arc_traits_2.h \
#    CGAL/Arr_circular_line_arc_traits_2.h \
#    CGAL/Arr_conic_traits_2.h \
#    CGAL/Arr_consolidated_curve_data_traits_2.h \
#    CGAL/Arr_counting_traits_2.h \
#    CGAL/Arr_curve_data_traits_2.h \
#    CGAL/Arr_dcel_base.h \
#    CGAL/Arr_default_dcel.h \
#    CGAL/Arr_default_overlay_traits.h \
#    CGAL/Arr_directional_non_caching_segment_basic_traits_2.h \
#    CGAL/Arr_enums.h \
#    CGAL/Arr_extended_dcel.h \
#    CGAL/Arr_face_index_map.h \
#    CGAL/Arr_face_map.h \
#    CGAL/Arr_geodesic_arc_on_sphere_partition_traits_2.h \
#    CGAL/Arr_geodesic_arc_on_sphere_traits_2.h \
#    CGAL/Arr_landmarks_point_location.h \
#    CGAL/Arr_line_arc_traits_2.h \
#    CGAL/Arr_linear_traits_2.h \
#    CGAL/Arr_naive_point_location.h \
#    CGAL/Arr_non_caching_segment_basic_traits_2.h \
#    CGAL/Arr_non_caching_segment_traits_2.h \
#    CGAL/Arr_observer.h \
#    CGAL/Arr_overlay.h \
#    CGAL/Arr_overlay_2.h \
#    CGAL/Arr_point_location_result.h \
#    CGAL/Arr_polycurve_basic_traits_2.h \
#    CGAL/Arr_polycurve_traits_2.h \
#    CGAL/Arr_polyline_traits_2.h \
#    CGAL/Arr_rational_function_traits_2.h \
#    CGAL/Arr_segment_traits_2.h \
#    CGAL/Arr_simple_point_location.h \
#    CGAL/Arr_spherical_topology_traits_2.h \
#    CGAL/Arr_tags.h \
#    CGAL/Arr_tracing_traits_2.h \
#    CGAL/Arr_trapezoid_ric_point_location.h \
#    CGAL/Arr_triangulation_point_location.h \
#    CGAL/Arr_unb_planar_topology_traits_2.h \
#    CGAL/Arr_vertex_index_map.h \
#    CGAL/Arr_vertex_map.h \
#    CGAL/Arr_vertical_decomposition_2.h \
#    CGAL/Arr_walk_along_line_point_location.h \
#    CGAL/arrange_offset_polygons_2.h \
#    CGAL/Arrangement_2.h \
#    CGAL/Arrangement_on_surface_2.h \
#    CGAL/Arrangement_on_surface_with_history_2.h \
#    CGAL/Arrangement_with_history_2.h \
#    CGAL/Arrangement_zone_2.h \
#    CGAL/array.h \
#    CGAL/assertions.h \
#    CGAL/assertions_behaviour.h \
#    CGAL/assertions_impl.h \
#    CGAL/atomic.h \
#    CGAL/barycenter.h \
#    CGAL/Barycentric_coordinates_2.h \
#    CGAL/basic.h \
#    CGAL/basic_classes.h \
#    CGAL/basic_constructions_2.h \
#    CGAL/basic_constructions_3.h \
    CGAL/Bbox_2.h \
    CGAL/Bbox_2_intersection.h \
    CGAL/Bbox_2_Line_2_intersection.h \
    CGAL/Bbox_2_Line_2_intersection_impl.h \
    CGAL/Bbox_2_Ray_2_intersection.h \
    CGAL/Bbox_3.h \
    CGAL/bbox_intersection_3.h \
#    CGAL/Bigfloat_interval_traits.h \
#    CGAL/bilateral_smooth_point_set.h \
#    CGAL/Boolean_set_operations_2.h \
#    CGAL/Bounded_kernel.h \
#    CGAL/bounding_box.h \
#    CGAL/box_intersection_d.h \
#    CGAL/Buffer_for_vao.h \
#    CGAL/Cache.h \
#    CGAL/Cartesian.h \
#    CGAL/Cartesian_converter.h \
#    CGAL/Cartesian_converter_fwd.h \
#    CGAL/Cartesian_d.h \
#    CGAL/cartesian_homogeneous_conversion.h \
#    CGAL/Cartesian_matrix.h \
#    CGAL/CC_safe_handle.h \
#    CGAL/Cell_attribute.h \
#    CGAL/Cell_attribute_with_id.h \
#    CGAL/Cell_attribute_with_point.h \
#    CGAL/Cell_attribute_with_point_and_id.h \
#    CGAL/Cell_const_iterators.h \
#    CGAL/Cell_iterators.h \
#    CGAL/centroid.h \
#    CGAL/certified_numeric_predicates.h \
#    CGAL/certified_quotient_predicates.h \
#    CGAL/CGAL_Ipelet_base.h \
#    CGAL/CGAL_Ipelet_base_v6.h \
#    CGAL/CGAL_Ipelet_base_v7.h \
#    CGAL/ch_akl_toussaint.h \
#    CGAL/ch_bykat.h \
#    CGAL/ch_eddy.h \
#    CGAL/ch_function_objects_2.h \
#    CGAL/ch_graham_andrew.h \
#    CGAL/ch_jarvis.h \
#    CGAL/ch_melkman.h \
#    CGAL/ch_selected_extreme_points_2.h \
#    CGAL/Chinese_remainder_traits.h \
#    CGAL/Circle_2.h \
#    CGAL/Circle_2_Circle_2_intersection.h \
#    CGAL/Circle_2_Line_2_intersection.h \
#    CGAL/Circle_3.h \
#    CGAL/Circle_type.h \
#    CGAL/Circular_arc_2.h \
#    CGAL/Circular_arc_3.h \
#    CGAL/Circular_arc_point_2.h \
#    CGAL/Circular_arc_point_3.h \
#    CGAL/Circular_kernel_2.h \
#    CGAL/Circular_kernel_converter.h \
#    CGAL/Circular_kernel_intersections.h \
#    CGAL/Circular_kernel_type_equality_wrapper.h \
#    CGAL/circulator.h \
#    CGAL/circulator_bases.h \
#    CGAL/Circulator_identity.h \
#    CGAL/Circulator_on_node.h \
#    CGAL/Circulator_project.h \
#    CGAL/Classification.h \
#    CGAL/CMap_linear_cell_complex_storages.h \
#    CGAL/Coercion_traits.h \
#    CGAL/Combination_enumerator.h \
#    CGAL/Combinatorial_map.h \
#    CGAL/Combinatorial_map_basic_operations.h \
#    CGAL/Combinatorial_map_constructors.h \
#    CGAL/Combinatorial_map_functors.h \
#    CGAL/Combinatorial_map_insertions.h \
#    CGAL/Combinatorial_map_iterators_base.h \
#    CGAL/Combinatorial_map_min_items.h \
#    CGAL/Combinatorial_map_operations.h \
#    CGAL/Combinatorial_map_save_load.h \
#    CGAL/Combinatorial_map_storages.h \
#    CGAL/Compact_container.h \
#    CGAL/Compact_mesh_cell_base_3.h \
#    CGAL/Compare_handles_with_or_without_timestamps.h \
#    CGAL/compare_vertices.h \
#    CGAL/Complex_2_in_triangulation_3.h \
#    CGAL/Complex_2_in_triangulation_cell_base_3.h \
#    CGAL/Complex_2_in_triangulation_vertex_base_3.h \
#    CGAL/Complexity_tags.h \
#    CGAL/Compute_anchor_3.h \
#    CGAL/compute_average_spacing.h \
#    CGAL/Compute_cone_boundaries_2.h \
#    CGAL/compute_outer_frame_margin.h \
#    CGAL/Concatenate_iterator.h \
#    CGAL/Concurrent_compact_container.h \
#    CGAL/Cone_spanners_enum_2.h \
#    CGAL/config.h \
#    CGAL/Conic_2.h \
#    CGAL/connect_holes.h \
#    CGAL/constant.h \
#    CGAL/Constrained_Delaunay_triangulation_2.h \
#    CGAL/Constrained_Delaunay_triangulation_face_base_2.h \
#    CGAL/Constrained_triangulation_2.h \
#    CGAL/Constrained_triangulation_face_base_2.h \
#    CGAL/Constrained_triangulation_plus_2.h \
#    CGAL/Constrained_voronoi_diagram_2.h \
#    CGAL/Constraint_hierarchy_2.h \
#    CGAL/Construct_theta_graph_2.h \
#    CGAL/Construct_yao_graph_2.h \
#    CGAL/constructions_d.h \
#    CGAL/convert_to_bfi.h \
#    CGAL/convex_decomposition_3.h \
#    CGAL/convex_hull_2.h \
#    CGAL/convex_hull_3.h \
#    CGAL/convex_hull_3_to_face_graph.h \
#    CGAL/convex_hull_3_to_polyhedron_3.h \
#    CGAL/convex_hull_constructive_traits_2.h \
#    CGAL/Convex_hull_d.h \
#    CGAL/Convex_hull_d_to_polyhedron_3.h \
#    CGAL/Convex_hull_d_traits_3.h \
#    CGAL/Convex_hull_face_base_2.h \
#    CGAL/Convex_hull_projective_xy_traits_2.h \
#    CGAL/Convex_hull_projective_xz_traits_2.h \
#    CGAL/Convex_hull_projective_yz_traits_2.h \
#    CGAL/convex_hull_traits_2.h \
#    CGAL/Convex_hull_traits_3.h \
#    CGAL/convexity_check_2.h \
#    CGAL/convexity_check_3.h \
#    CGAL/copy_n.h \
#    CGAL/CORE_algebraic_number_traits.h \
#    CGAL/CORE_arithmetic_kernel.h \
#    CGAL/CORE_BigFloat.h \
#    CGAL/CORE_BigInt.h \
#    CGAL/CORE_BigRat.h \
#    CGAL/CORE_coercion_traits.h \
#    CGAL/CORE_Expr.h \
#    CGAL/corefinement_operations.h \
#    CGAL/Counted_number.h \
#    CGAL/Counting_iterator.h \
#    CGAL/create_offset_polygons_2.h \
#    CGAL/create_offset_polygons_from_polygon_with_holes_2.h \
#    CGAL/create_straight_skeleton_2.h \
#    CGAL/create_straight_skeleton_from_polygon_with_holes_2.h \
#    CGAL/Dart.h \
#    CGAL/Dart_const_iterators.h \
#    CGAL/Dart_iterators.h \
#    CGAL/Default.h \
#    CGAL/Default_diagonalize_traits.h \
#    CGAL/Deformation_Eigen_closest_rotation_traits_3.h \
#    CGAL/Deformation_Eigen_polar_closest_rotation_traits_3.h \
#    CGAL/Delaunay_d.h \
#    CGAL/Delaunay_mesh_area_criteria_2.h \
#    CGAL/Delaunay_mesh_criteria_2.h \
#    CGAL/Delaunay_mesh_face_base_2.h \
#    CGAL/Delaunay_mesh_local_size_criteria_2.h \
#    CGAL/Delaunay_mesh_size_criteria_2.h \
#    CGAL/Delaunay_mesh_vertex_base_2.h \
#    CGAL/Delaunay_mesher_2.h \
#    CGAL/Delaunay_mesher_no_edge_refinement_2.h \
#    CGAL/Delaunay_triangulation.h \
#    CGAL/Delaunay_triangulation_2.h \
#    CGAL/Delaunay_triangulation_3.h \
#    CGAL/Delaunay_triangulation_adaptation_policies_2.h \
#    CGAL/Delaunay_triangulation_adaptation_traits_2.h \
#    CGAL/Delaunay_triangulation_cell_base_3.h \
#    CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h \
#    CGAL/demangle.h \
#    CGAL/determinant.h \
#    CGAL/determinant_of_vectors.h \
#    CGAL/Diagonalize_traits.h \
#    CGAL/Dimension.h \
#    CGAL/Direction_2.h \
#    CGAL/Direction_3.h \
#    CGAL/disable_warnings.h \
#    CGAL/Distance_2.h \
#    CGAL/distance_predicates_2.h \
#    CGAL/distance_predicates_3.h \
#    CGAL/double.h \
#    CGAL/Double_map.h \
#    CGAL/draw_linear_cell_complex.h \
#    CGAL/draw_polyhedron.h \
#    CGAL/draw_surface_mesh.h \
#    CGAL/draw_triangulation_2.h \
#    CGAL/draw_triangulation_3.h \
#    CGAL/Dummy_tds_2.h \
#    CGAL/Dynamic_matrix.h \
#    CGAL/Dynamic_property_map.h \
#    CGAL/edge_aware_upsample_point_set.h \
#    CGAL/Eigen_diagonalize_traits.h \
#    CGAL/Eigen_matrix.h \
#    CGAL/Eigen_solver_traits.h \
#    CGAL/Eigen_svd.h \
#    CGAL/Eigen_vector.h \
#    CGAL/enable_warnings.h \
#    CGAL/enum.h \
#    CGAL/Enum_converter.h \
#    CGAL/Env_default_diagram_1.h \
#    CGAL/Env_plane_traits_3.h \
#    CGAL/Env_sphere_traits_3.h \
#    CGAL/Env_surface_data_traits_3.h \
#    CGAL/Env_tracing_traits_3.h \
#    CGAL/Env_triangle_traits_3.h \
#    CGAL/envelope_2.h \
#    CGAL/envelope_3.h \
#    CGAL/Envelope_diagram_1.h \
#    CGAL/Epeck_d.h \
#    CGAL/Epic_converter.h \
#    CGAL/Epick_d.h \
#    CGAL/estimate_scale.h \
#    CGAL/Euclidean_distance.h \
#    CGAL/Euclidean_distance_sphere_point.h \
#    CGAL/Euler_integrator_2.h \
#    CGAL/Exact_circular_kernel_2.h \
#    CGAL/Exact_integer.h \
#    CGAL/Exact_kernel_selector.h \
#    CGAL/Exact_predicates_exact_constructions_kernel.h \
#    CGAL/Exact_predicates_exact_constructions_kernel_with_kth_root.h \
#    CGAL/Exact_predicates_exact_constructions_kernel_with_root_of.h \
#    CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h \
#    CGAL/Exact_predicates_inexact_constructions_kernel.h \
#    CGAL/Exact_rational.h \
#    CGAL/Exact_spherical_kernel_3.h \
#    CGAL/exceptions.h \
#    CGAL/Exponent_vector.h \
#    CGAL/Extended_cartesian.h \
#    CGAL/extended_euclidean_algorithm.h \
#    CGAL/Extended_homogeneous.h \
#    CGAL/extract_mean_curvature_flow_skeleton.h \
#    CGAL/extremal_polygon_2.h \
#    CGAL/Extremal_polygon_traits_2.h \
#    CGAL/Extreme_points_traits_adapter_3.h \
#    CGAL/exude_mesh_3.h \
#    CGAL/Filter_circulator.h \
#    CGAL/Filtered_bbox_circular_kernel_2.h \
#    CGAL/Filtered_construction.h \
#    CGAL/Filtered_extended_homogeneous.h \
#    CGAL/Filtered_kernel.h \
#    CGAL/Filtered_kernel_d.h \
#    CGAL/Filtered_kernel_fwd.h \
#    CGAL/Filtered_predicate.h \
#    CGAL/Filtered_predicate_with_state.h \
#    CGAL/Fixed_alpha_shape_3.h \
#    CGAL/Fixed_alpha_shape_cell_base_3.h \
#    CGAL/Fixed_alpha_shape_vertex_base_3.h \
#    CGAL/Flattening_iterator.h \
#    CGAL/float.h \
#    CGAL/Fourtuple.h \
#    CGAL/FPU.h \
#    CGAL/FPU_extension.h \
#    CGAL/FPU_gcc_i386.h \
#    CGAL/FPU_gcc_i386_sse2.h \
#    CGAL/FPU_msvc.h \
#    CGAL/Fraction_traits.h \
#    CGAL/function.h \
#    CGAL/function_objects.h \
#    CGAL/functional.h \
#    CGAL/functions_on_enums.h \
#    CGAL/functions_on_signs.h \
#    CGAL/Fuzzy_iso_box.h \
#    CGAL/Fuzzy_sphere.h \
#    CGAL/General_polygon_2.h \
#    CGAL/General_polygon_set_2.h \
#    CGAL/General_polygon_set_on_surface_2.h \
#    CGAL/General_polygon_with_holes_2.h \
#    CGAL/Generalized_map.h \
#    CGAL/Generalized_map_iterators_base.h \
#    CGAL/Generalized_map_operations.h \
#    CGAL/Generalized_map_save_load.h \
#    CGAL/Generalized_map_storages.h \
#    CGAL/generators.h \
#    CGAL/Generic_map_min_items.h \
#    CGAL/generic_sweep.h \
#    CGAL/Get_arithmetic_kernel.h \
#    CGAL/gl.h \
#    CGAL/global_functions_circular_kernel_2.h \
#    CGAL/global_functions_on_root_for_sphere_2_3.h \
#    CGAL/global_functions_on_roots_and_polynomials_1_3.h \
#    CGAL/global_functions_on_roots_and_polynomials_2_3.h \
#    CGAL/global_functions_spherical_kernel_3.h \
#    CGAL/glu.h \
#    CGAL/GMap_cell_const_iterators.h \
#    CGAL/GMap_cell_iterators.h \
#    CGAL/GMap_dart_const_iterators.h \
#    CGAL/GMap_dart_iterators.h \
#    CGAL/GMap_linear_cell_complex_storages.h \
#    CGAL/gmp.h \
#    CGAL/GMP_arithmetic_kernel.h \
#    CGAL/Gmp_coercion_traits.h \
#    CGAL/Gmpfi.h \
#    CGAL/Gmpfr.h \
#    CGAL/Gmpq.h \
#    CGAL/gmpxx.h \
#    CGAL/GMPXX_arithmetic_kernel.h \
#    CGAL/gmpxx_coercion_traits.h \
#    CGAL/Gmpz.h \
#    CGAL/Gmpzf.h \
#    CGAL/gnuplot_output_2.h \
#    CGAL/Gps_circle_segment_traits_2.h \
#    CGAL/Gps_segment_traits_2.h \
#    CGAL/Gps_traits_2.h \
#    CGAL/grabbers.h \
#    CGAL/graph_traits_Arrangement_2.h \
#    CGAL/graph_traits_dual_arrangement_2.h \
#    CGAL/graph_traits_dual_arrangement_on_surface_2.h \
#    CGAL/graph_traits_dual_arrangement_on_surface_with_history_2.h \
#    CGAL/graph_traits_dual_arrangement_with_history_2.h \
#    CGAL/Gray_image_mesh_domain_3.h \
#    CGAL/Gray_level_image_3.h \
#    CGAL/grid_simplify_point_set.h \
#    CGAL/halfedgeds_connected_components.h \
#    CGAL/HalfedgeDS_const_decorator.h \
#    CGAL/halfedgeDS_cut_component.h \
#    CGAL/HalfedgeDS_decorator.h \
#    CGAL/HalfedgeDS_default.h \
#    CGAL/HalfedgeDS_face_base.h \
#    CGAL/HalfedgeDS_face_max_base_with_id.h \
#    CGAL/HalfedgeDS_face_min_base.h \
#    CGAL/HalfedgeDS_halfedge_base.h \
#    CGAL/HalfedgeDS_halfedge_max_base_with_id.h \
#    CGAL/HalfedgeDS_halfedge_min_base.h \
#    CGAL/HalfedgeDS_items_2.h \
#    CGAL/HalfedgeDS_items_decorator.h \
#    CGAL/HalfedgeDS_iterator.h \
#    CGAL/HalfedgeDS_iterator_adaptor.h \
#    CGAL/HalfedgeDS_list.h \
#    CGAL/HalfedgeDS_min_items.h \
#    CGAL/HalfedgeDS_vector.h \
#    CGAL/HalfedgeDS_vertex_base.h \
#    CGAL/HalfedgeDS_vertex_max_base_with_id.h \
#    CGAL/HalfedgeDS_vertex_min_base.h \
#    CGAL/Handle.h \
#    CGAL/Handle_for.h \
#    CGAL/Handle_for_virtual.h \
#    CGAL/Handle_hash_function.h \
#    CGAL/Handle_with_policy.h \
#    CGAL/Has_conversion.h \
#    CGAL/Has_timestamp.h \
#    CGAL/Hash_handles_with_or_without_timestamps.h \
#    CGAL/hash_openmesh.h \
#    CGAL/Hidden_point_memory_policy.h \
#    CGAL/hierarchy_simplify_point_set.h \
#    CGAL/Hilbert_policy_tags.h \
#    CGAL/hilbert_sort.h \
#    CGAL/Hilbert_sort_2.h \
#    CGAL/Hilbert_sort_3.h \
#    CGAL/Hilbert_sort_base.h \
#    CGAL/Hilbert_sort_d.h \
#    CGAL/Hilbert_sort_median_2.h \
#    CGAL/Hilbert_sort_median_3.h \
#    CGAL/Hilbert_sort_median_d.h \
#    CGAL/Hilbert_sort_middle_2.h \
#    CGAL/Hilbert_sort_middle_3.h \
#    CGAL/Hilbert_sort_middle_base.h \
#    CGAL/Hilbert_sort_middle_d.h \
#    CGAL/hilbert_sort_on_sphere.h \
#    CGAL/Hilbert_sort_on_sphere_3.h \
#    CGAL/Homogeneous.h \
#    CGAL/Homogeneous_converter.h \
#    CGAL/Homogeneous_d.h \
#    CGAL/Hyperbola_2.h \
#    CGAL/Hyperbola_ray_2.h \
#    CGAL/Hyperbola_segment_2.h \
#    CGAL/Identity_policy_2.h \
#    CGAL/IEEE_754_unions.h \
#    CGAL/Image_3.h \
#    CGAL/Image_3_impl.h \
#    CGAL/Image_3_vtk_interface.h \
#    CGAL/ImageIO.h \
#    CGAL/ImageIO_impl.h \
#    CGAL/Implicit_mesh_domain_3.h \
#    CGAL/Implicit_surface_3.h \
#    CGAL/Implicit_to_labeled_subdomains_function_wrapper.h \
#    CGAL/Implicit_to_labeling_function_wrapper.h \
#    CGAL/in_place_edge_list.h \
#    CGAL/In_place_list.h \
#    CGAL/Incremental_neighbor_search.h \
#    CGAL/Index_property_map.h \
#    CGAL/int.h \
#    CGAL/interpolation_functions.h \
#    CGAL/Interpolation_gradient_fitting_traits_2.h \
#    CGAL/Interpolation_traits_2.h \
#    CGAL/intersection_2.h \
#    CGAL/intersection_2_1.h \
#    CGAL/intersection_2_2.h \
#    CGAL/intersection_2_3.h \
#    CGAL/intersection_3.h \
#    CGAL/intersection_3_1.h \
#    CGAL/intersection_of_Polyhedra_3.h \
#    CGAL/intersection_of_Polyhedra_3_refinement_visitor.h \
#    CGAL/Intersection_traits.h \
#    CGAL/Intersection_traits_2.h \
#    CGAL/Intersection_traits_3.h \
#    CGAL/intersections.h \
#    CGAL/intersections_d.h \
#    CGAL/Interval_arithmetic.h \
#    CGAL/Interval_arithmetic_impl.h \
#    CGAL/Interval_nt.h \
#    CGAL/Interval_skip_list.h \
#    CGAL/Interval_skip_list_interval.h \
#    CGAL/Interval_traits.h \
#    CGAL/Inverse_index.h \
#    CGAL/ipower.h \
#    CGAL/Is_a_predicate.h \
#    CGAL/is_convertible.h \
#    CGAL/Is_extended_kernel.h \
#    CGAL/is_iterator.h \
#    CGAL/is_streamable.h \
#    CGAL/is_y_monotone_2.h \
#    CGAL/Iso_cuboid_3.h \
#    CGAL/Iso_rectangle_2.h \
#    CGAL/Iso_rectangle_2_Iso_rectangle_2_intersection.h \
#    CGAL/Iso_rectangle_2_Line_2_intersection.h \
#    CGAL/Iso_rectangle_2_Point_2_intersection.h \
#    CGAL/Iso_rectangle_2_Ray_2_intersection.h \
#    CGAL/Iso_rectangle_2_Segment_2_intersection.h \
#    CGAL/Iso_rectangle_d.h \
#    CGAL/iterator.h \
#    CGAL/iterator_from_indices.h \
#    CGAL/Iterator_project.h \
#    CGAL/Iterator_range.h \
#    CGAL/Iterator_transform.h \
#    CGAL/jet_estimate_normals.h \
#    CGAL/jet_smooth_point_set.h \
#    CGAL/Join_input_iterator.h \
#    CGAL/K_neighbor_search.h \
#    CGAL/Kd_tree.h \
#    CGAL/Kd_tree_node.h \
#    CGAL/Kd_tree_rectangle.h \
#    CGAL/kernel_assertions.h \
#    CGAL/kernel_basic.h \
#    CGAL/Kernel_checker.h \
#    CGAL/kernel_config.h \
#    CGAL/Kernel_profiler.h \
#    CGAL/kernel_to_kernel.h \
#    CGAL/Kernel_traits.h \
#    CGAL/Kernel_traits_fwd.h \
#    CGAL/known_bit_size_integers.h \
#    CGAL/Labeled_image_mesh_domain_3.h \
#    CGAL/Labeled_mesh_domain_3.h \
#    CGAL/Lapack_svd.h \
#    CGAL/Largest_empty_iso_rectangle_2.h \
#    CGAL/Lazy.h \
#    CGAL/Lazy_exact_nt.h \
#    CGAL/Lazy_kernel.h \
#    CGAL/LEDA_arithmetic_kernel.h \
#    CGAL/LEDA_basic.h \
#    CGAL/leda_bigfloat.h \
#    CGAL/leda_bigfloat_interval.h \
#    CGAL/leda_coercion_traits.h \
#    CGAL/leda_integer.h \
#    CGAL/leda_rational.h \
#    CGAL/leda_real.h \
#    CGAL/Level_interval.h \
#    CGAL/license.h \
#    CGAL/Lightweight_vector_3.h \
#    CGAL/Line_2.h \
#    CGAL/Line_2_Bbox_2_intersection.h \
#    CGAL/Line_2_Iso_rectangle_2_intersection.h \
#    CGAL/Line_2_Line_2_intersection.h \
#    CGAL/Line_2_Point_2_intersection.h \
#    CGAL/Line_2_Ray_2_intersection.h \
#    CGAL/Line_2_Segment_2_intersection.h \
#    CGAL/Line_2_Triangle_2_intersection.h \
#    CGAL/Line_3.h \
#    CGAL/Line_arc_2.h \
#    CGAL/Line_arc_3.h \
#    CGAL/Linear_algebraCd.h \
#    CGAL/Linear_algebraHd.h \
#    CGAL/Linear_cell_complex.h \
#    CGAL/Linear_cell_complex_base.h \
#    CGAL/Linear_cell_complex_bgl_min_items.h \
#    CGAL/Linear_cell_complex_constructors.h \
#    CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h \
#    CGAL/Linear_cell_complex_for_combinatorial_map.h \
#    CGAL/Linear_cell_complex_for_generalized_map.h \
#    CGAL/Linear_cell_complex_incremental_builder.h \
#    CGAL/Linear_cell_complex_min_items.h \
#    CGAL/Linear_cell_complex_operations.h \
#    CGAL/Linear_cell_complex_traits.h \
#    CGAL/linear_least_squares_fitting_2.h \
#    CGAL/linear_least_squares_fitting_3.h \
#    CGAL/linear_least_squares_fitting_circles_2.h \
#    CGAL/linear_least_squares_fitting_cuboids_3.h \
#    CGAL/linear_least_squares_fitting_points_2.h \
#    CGAL/linear_least_squares_fitting_points_3.h \
#    CGAL/linear_least_squares_fitting_rectangles_2.h \
#    CGAL/linear_least_squares_fitting_segments_2.h \
#    CGAL/linear_least_squares_fitting_segments_3.h \
#    CGAL/linear_least_squares_fitting_spheres_3.h \
#    CGAL/linear_least_squares_fitting_tetrahedra_3.h \
#    CGAL/linear_least_squares_fitting_triangles_2.h \
#    CGAL/linear_least_squares_fitting_triangles_3.h \
#    CGAL/link_to_face_graph.h \
#    CGAL/lloyd_optimize_mesh_2.h \
#    CGAL/lloyd_optimize_mesh_3.h \
#    CGAL/Location_policy.h \
#    CGAL/long_double.h \
#    CGAL/long_long.h \
#    CGAL/make_mesh_3.h \
#    CGAL/make_periodic_3_mesh_3.h \
#    CGAL/make_piecewise_smooth_surface_mesh.h \
#    CGAL/make_skin_surface_mesh_3.h \
#    CGAL/make_surface_mesh.h \
#    CGAL/make_union_of_balls_3.h \
#    CGAL/Manhattan_distance_iso_box_point.h \
#    CGAL/marching_tetrahedra_3.h \
#    CGAL/Marching_tetrahedra_observer_default_3.h \
#    CGAL/Marching_tetrahedra_traits_skin_surface_3.h \
#    CGAL/Mean_curvature_flow_skeletonization.h \
#    CGAL/memory.h \
#    CGAL/Memory_sizer.h \
#    CGAL/Mesh_cell_base_3.h \
#    CGAL/Mesh_cell_criteria_3.h \
#    CGAL/Mesh_complex_3_in_triangulation_3.h \
#    CGAL/Mesh_constant_domain_field_3.h \
#    CGAL/Mesh_criteria_3.h \
#    CGAL/Mesh_domain_with_polyline_features_3.h \
#    CGAL/Mesh_edge_criteria_3.h \
#    CGAL/Mesh_error_code.h \
#    CGAL/Mesh_facet_criteria_3.h \
#    CGAL/Mesh_facet_topology.h \
#    CGAL/Mesh_optimization_return_code.h \
#    CGAL/Mesh_polyhedron_3.h \
#    CGAL/mesh_segmentation.h \
#    CGAL/mesh_skin_surface_3.h \
#    CGAL/Mesh_triangulation_3.h \
#    CGAL/mesh_union_of_balls_3.h \
#    CGAL/Mesh_vertex_base_3.h \
#    CGAL/Mesher_level.h \
#    CGAL/Mesher_level_default_implementations.h \
#    CGAL/Mesher_level_visitors.h \
#    CGAL/Min_annulus_d.h \
#    CGAL/Min_circle_2.h \
#    CGAL/Min_circle_2_traits_2.h \
#    CGAL/Min_ellipse_2.h \
#    CGAL/Min_ellipse_2_traits_2.h \
#    CGAL/min_max_n.h \
#    CGAL/min_quadrilateral_2.h \
#    CGAL/Min_quadrilateral_traits_2.h \
#    CGAL/Min_sphere_annulus_d_traits_2.h \
#    CGAL/Min_sphere_annulus_d_traits_3.h \
#    CGAL/Min_sphere_annulus_d_traits_d.h \
#    CGAL/Min_sphere_d.h \
#    CGAL/Min_sphere_of_points_d_traits_2.h \
#    CGAL/Min_sphere_of_points_d_traits_3.h \
#    CGAL/Min_sphere_of_points_d_traits_d.h \
#    CGAL/Min_sphere_of_spheres_d.h \
#    CGAL/Min_sphere_of_spheres_d_traits_2.h \
#    CGAL/Min_sphere_of_spheres_d_traits_3.h \
#    CGAL/Min_sphere_of_spheres_d_traits_d.h \
#    CGAL/minimum_enclosing_quadrilateral_2.h \
#    CGAL/Minimum_enclosing_quadrilateral_traits_2.h \
#    CGAL/minkowski_sum_2.h \
#    CGAL/minkowski_sum_3.h \
#    CGAL/Modifiable_priority_queue.h \
#    CGAL/Modifier_base.h \
#    CGAL/Modular_traits.h \
#    CGAL/Monge_via_jet_fitting.h \
#    CGAL/monotone_matrix_search.h \
#    CGAL/more_functions_on_signs.h \
#    CGAL/MP_Float.h \
#    CGAL/MP_Float_arithmetic_kernel.h \
#    CGAL/MP_Float_impl.h \
#    CGAL/mpfi_coercion_traits.h \
#    CGAL/mpfr_coercion_traits.h \
#    CGAL/mpq_class.h \
#    CGAL/mpz_class.h \
#    CGAL/Mpzf.h \
#    CGAL/mst_orient_normals.h \
#    CGAL/MSVC_compiler_config.h \
#    CGAL/Multi_surface_3.h \
#    CGAL/Multiscale_sort.h \
#    CGAL/Multiset.h \
#    CGAL/multiset_assertions.h \
#    CGAL/mutex.h \
#    CGAL/N_step_adaptor.h \
#    CGAL/N_step_adaptor_derived.h \
#    CGAL/natural_neighbor_coordinates_2.h \
#    CGAL/natural_neighbor_coordinates_3.h \
#    CGAL/nearest_neighbor_delaunay_2.h \
#    CGAL/Needs_parens_as_product.h \
#    CGAL/Nef_nary_intersection_3.h \
#    CGAL/Nef_nary_union_3.h \
#    CGAL/Nef_polyhedron_2.h \
#    CGAL/Nef_polyhedron_3.h \
#    CGAL/Nef_polyhedron_S2.h \
#    CGAL/Nef_polynomial.h \
#    CGAL/Nef_polynomial_fwd.h \
#    CGAL/Nested_iterator.h \
#    CGAL/No_intersection_surface_sweep_2.h \
#    CGAL/normal_vector_newell_3.h \
#    CGAL/NT_converter.h \
#    CGAL/Null_matrix.h \
#    CGAL/number_type_basic.h \
#    CGAL/Number_type_checker.h \
#    CGAL/number_type_config.h \
#    CGAL/number_utils.h \
#    CGAL/number_utils_classes.h \
#    CGAL/Object.h \
#    CGAL/odt_optimize_mesh_3.h \
#    CGAL/OFF_to_nef_3.h \
#    CGAL/offset_polygon_2.h \
#    CGAL/Optimal_transportation_reconstruction_2.h \
#    CGAL/Optimisation_d_traits_2.h \
#    CGAL/Optimisation_d_traits_3.h \
#    CGAL/Optimisation_d_traits_d.h \
#    CGAL/optimize_mesh_3.h \
#    CGAL/optimize_periodic_3_mesh_3.h \
#    CGAL/Orientation_Linf_2.h \
#    CGAL/Origin.h \
#    CGAL/Origin_impl.h \
#    CGAL/Orthogonal_incremental_neighbor_search.h \
#    CGAL/Orthogonal_k_neighbor_search.h \
#    CGAL/Parabola_2.h \
#    CGAL/Parabola_segment_2.h \
#    CGAL/partition_2.h \
#    CGAL/partition_is_valid_2.h \
#    CGAL/Partition_is_valid_traits_2.h \
#    CGAL/Partition_traits_2.h \
#    CGAL/pca_estimate_normals.h \
#    CGAL/PCA_util.h \
#    CGAL/PCA_util_Eigen.h \
#    CGAL/Periodic_2_Delaunay_triangulation_2.h \
#    CGAL/Periodic_2_Delaunay_triangulation_traits_2.h \
#    CGAL/Periodic_2_offset_2.h \
#    CGAL/Periodic_2_triangulation_2.h \
#    CGAL/Periodic_2_triangulation_dummy_12.h \
#    CGAL/Periodic_2_triangulation_face_base_2.h \
#    CGAL/Periodic_2_triangulation_hierarchy_2.h \
#    CGAL/Periodic_2_triangulation_iterators_2.h \
#    CGAL/Periodic_2_triangulation_traits_2.h \
#    CGAL/Periodic_2_triangulation_vertex_base_2.h \
#    CGAL/Periodic_3_Delaunay_triangulation_3.h \
#    CGAL/Periodic_3_Delaunay_triangulation_traits_3.h \
#    CGAL/Periodic_3_function_wrapper.h \
#    CGAL/Periodic_3_mesh_triangulation_3.h \
#    CGAL/Periodic_3_offset_3.h \
#    CGAL/Periodic_3_regular_triangulation_3.h \
#    CGAL/Periodic_3_regular_triangulation_traits_3.h \
#    CGAL/Periodic_3_triangulation_3.h \
#    CGAL/periodic_3_triangulation_3_io.h \
#    CGAL/Periodic_3_triangulation_ds_cell_base_3.h \
#    CGAL/Periodic_3_triangulation_ds_vertex_base_3.h \
#    CGAL/Periodic_3_triangulation_hierarchy_3.h \
#    CGAL/Periodic_3_triangulation_traits_3.h \
#    CGAL/perturb_mesh_3.h \
#    CGAL/pierce_rectangles_2.h \
#    CGAL/Plane_3.h \
#    CGAL/Plane_separator.h \
#    CGAL/Point_2.h \
#    CGAL/Point_2_Iso_rectangle_2_intersection.h \
#    CGAL/Point_2_Line_2_intersection.h \
#    CGAL/Point_2_Point_2_intersection.h \
#    CGAL/Point_2_Ray_2_intersection.h \
#    CGAL/Point_2_Segment_2_intersection.h \
#    CGAL/Point_2_Triangle_2_intersection.h \
#    CGAL/Point_3.h \
#    CGAL/Point_container.h \
#    CGAL/point_generators_2.h \
#    CGAL/point_generators_3.h \
#    CGAL/point_generators_d.h \
#    CGAL/Point_set_2.h \
#    CGAL/Point_set_3.h \
#    CGAL/point_set_processing_assertions.h \
#    CGAL/Point_traits.h \
#    CGAL/Point_with_normal_3.h \
#    CGAL/Point_with_psc_localisation.h \
#    CGAL/Point_with_surface_index.h \
#    CGAL/Point_with_surface_index_geom_traits.h \
#    CGAL/Poisson_implicit_surface_3.h \
#    CGAL/Poisson_mesh_cell_criteria_3.h \
#    CGAL/Poisson_reconstruction_function.h \
#    CGAL/poisson_refine_triangulation.h \
#    CGAL/poisson_surface_reconstruction.h \
#    CGAL/Polychain_2.h \
#    CGAL/Polygon_2.h \
#    CGAL/Polygon_2_algorithms.h \
#    CGAL/Polygon_convex_decomposition_2.h \
#    CGAL/polygon_function_objects.h \
#    CGAL/polygon_mesh_processing.h \
#    CGAL/Polygon_mesh_slicer.h \
#    CGAL/Polygon_nop_decomposition_2.h \
#    CGAL/Polygon_offset_builder_2.h \
#    CGAL/Polygon_offset_builder_traits_2.h \
#    CGAL/Polygon_set_2.h \
#    CGAL/Polygon_traits_2.h \
#    CGAL/Polygon_triangulation_decomposition_2.h \
#    CGAL/Polygon_vertical_decomposition_2.h \
#    CGAL/Polygon_with_holes_2.h \
#    CGAL/Polyhedral_complex_mesh_domain_3.h \
#    CGAL/Polyhedral_mesh_domain_3.h \
#    CGAL/Polyhedral_mesh_domain_with_features_3.h \
#    CGAL/PolyhedralSurf_neighbors.h \
#    CGAL/Polyhedron_3.h \
#    CGAL/Polyhedron_3_fwd.h \
#    CGAL/Polyhedron_3_to_lcc.h \
#    CGAL/Polyhedron_copy_3.h \
#    CGAL/polyhedron_cut_plane_3.h \
#    CGAL/Polyhedron_incremental_builder_3.h \
#    CGAL/Polyhedron_items_3.h \
#    CGAL/Polyhedron_items_with_id_3.h \
#    CGAL/Polyhedron_min_items_3.h \
#    CGAL/Polyhedron_traits_3.h \
#    CGAL/Polyhedron_traits_with_normals_3.h \
#    CGAL/Polyline_constraint_hierarchy_2.h \
#    CGAL/Polynomial.h \
#    CGAL/Polynomial_traits_d.h \
#    CGAL/Polynomial_type_generator.h \
#    CGAL/polynomial_utils.h \
#    CGAL/Polynomials_1_2.h \
#    CGAL/Polynomials_1_3.h \
#    CGAL/Polynomials_2_2.h \
#    CGAL/Polynomials_2_3.h \
#    CGAL/Polynomials_for_line_3.h \
#    CGAL/Polytope_distance_d.h \
#    CGAL/Polytope_distance_d_traits_2.h \
#    CGAL/Polytope_distance_d_traits_3.h \
#    CGAL/Polytope_distance_d_traits_d.h \
#    CGAL/predicates_d.h \
#    CGAL/predicates_on_lines_2.h \
#    CGAL/predicates_on_points_2.h \
#    CGAL/predicates_on_points_3.h \
#    CGAL/primes.h \
#    CGAL/Profile_counter.h \
#    CGAL/Profile_timer.h \
#    CGAL/Projection_traits_xy_3.h \
#    CGAL/Projection_traits_xz_3.h \
#    CGAL/Projection_traits_yz_3.h \
#    CGAL/property_map.h \
#    CGAL/QP_functions.h \
#    CGAL/QP_models.h \
#    CGAL/QP_options.h \
#    CGAL/QP_solution.h \
#    CGAL/Quotient.h \
#    CGAL/Quotient_fwd.h \
#    CGAL/radial_orient_normals.h \
#    CGAL/Random.h \
#    CGAL/Random_access_adaptor.h \
#    CGAL/Random_access_value_adaptor.h \
#    CGAL/random_convex_hull_in_disc_2.h \
#    CGAL/Random_convex_hull_traits_2.h \
#    CGAL/random_convex_set_2.h \
#    CGAL/Random_convex_set_traits_2.h \
#    CGAL/Random_impl.h \
#    CGAL/random_polygon_2.h \
#    CGAL/Random_polygon_2_sweep.h \
#    CGAL/Random_polygon_traits_2.h \
#    CGAL/random_selection.h \
#    CGAL/random_simplify_point_set.h \
#    CGAL/range_search_delaunay_2.h \
#    CGAL/Range_segment_tree_traits.h \
#    CGAL/Range_tree_d.h \
#    CGAL/Range_tree_k.h \
#    CGAL/rational_rotation.h \
#    CGAL/Rational_traits.h \
#    CGAL/Ray_2.h \
#    CGAL/Ray_2_Bbox_2_intersection.h \
#    CGAL/Ray_2_Bbox_2_intersection_impl.h \
#    CGAL/Ray_2_Iso_rectangle_2_intersection.h \
#    CGAL/Ray_2_Line_2_intersection.h \
#    CGAL/Ray_2_Point_2_intersection.h \
#    CGAL/Ray_2_Ray_2_intersection.h \
#    CGAL/Ray_2_Segment_2_intersection.h \
#    CGAL/Ray_2_Triangle_2_intersection.h \
#    CGAL/Ray_3.h \
#    CGAL/read_vtk_image_data.h \
#    CGAL/Real_embeddable_traits.h \
#    CGAL/Real_timer.h \
#    CGAL/Real_timer_impl.h \
#    CGAL/Reconstruction_triangulation_3.h \
#    CGAL/rectangular_3_center_2.h \
#    CGAL/rectangular_p_center_2.h \
#    CGAL/Rectangular_p_center_traits_2.h \
#    CGAL/Referenced_argument.h \
#    CGAL/refine_mesh_3.h \
#    CGAL/refine_periodic_3_mesh_3.h \
#    CGAL/Regular_complex_d.h \
#    CGAL/Regular_grid_2.h \
#    CGAL/regular_neighbor_coordinates_2.h \
#    CGAL/Regular_triangulation.h \
#    CGAL/Regular_triangulation_2.h \
#    CGAL/Regular_triangulation_3.h \
#    CGAL/Regular_triangulation_adaptation_policies_2.h \
#    CGAL/Regular_triangulation_adaptation_traits_2.h \
#    CGAL/Regular_triangulation_cell_base_3.h \
#    CGAL/Regular_triangulation_cell_base_with_weighted_circumcenter_3.h \
#    CGAL/Regular_triangulation_euclidean_traits_2.h \
#    CGAL/Regular_triangulation_euclidean_traits_3.h \
#    CGAL/Regular_triangulation_face_base_2.h \
#    CGAL/Regular_triangulation_filtered_traits_2.h \
#    CGAL/Regular_triangulation_traits_adapter.h \
#    CGAL/Regular_triangulation_vertex_base_2.h \
#    CGAL/Regular_triangulation_vertex_base_3.h \
#    CGAL/regularize_planes.h \
#    CGAL/remove_far_points_in_mesh_3.h \
#    CGAL/remove_outliers.h \
#    CGAL/representation_tags.h \
#    CGAL/Residue.h \
#    CGAL/result_of.h \
#    CGAL/Ridges.h \
#    CGAL/Robust_circumcenter_filtered_traits_3.h \
#    CGAL/Robust_circumcenter_traits_3.h \
#    CGAL/Robust_construction.h \
#    CGAL/Robust_weighted_circumcenter_filtered_traits_3.h \
#    CGAL/Root_for_circles_2_2.h \
#    CGAL/Root_for_spheres_2_3.h \
#    CGAL/Root_of_traits.h \
#    CGAL/Root_of_traits_specializations.h \
#    CGAL/Rotational_sweep_visibility_2.h \
#    CGAL/Runge_kutta_integrator_2.h \
#    CGAL/Scalar_factor_traits.h \
#    CGAL/Scale_space_surface_reconstruction_3.h \
#    CGAL/Search_traits.h \
#    CGAL/Search_traits_2.h \
#    CGAL/Search_traits_3.h \
#    CGAL/Search_traits_adapter.h \
#    CGAL/Search_traits_d.h \
#    CGAL/Search_traits_vertex_handle_3.h \
#    CGAL/Segment_2.h \
#    CGAL/Segment_2_Iso_rectangle_2_intersection.h \
#    CGAL/Segment_2_Line_2_intersection.h \
#    CGAL/Segment_2_Point_2_intersection.h \
#    CGAL/Segment_2_Ray_2_intersection.h \
#    CGAL/Segment_2_Segment_2_intersection.h \
#    CGAL/Segment_2_Triangle_2_intersection.h \
#    CGAL/Segment_3.h \
#    CGAL/Segment_Delaunay_graph_2.h \
#    CGAL/Segment_Delaunay_graph_adaptation_policies_2.h \
#    CGAL/Segment_Delaunay_graph_adaptation_traits_2.h \
#    CGAL/Segment_Delaunay_graph_face_base_2.h \
#    CGAL/Segment_Delaunay_graph_filtered_traits_2.h \
#    CGAL/Segment_Delaunay_graph_hierarchy_2.h \
#    CGAL/Segment_Delaunay_graph_hierarchy_vertex_base_2.h \
#    CGAL/Segment_Delaunay_graph_Linf_2.h \
#    CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h \
#    CGAL/Segment_Delaunay_graph_Linf_hierarchy_2.h \
#    CGAL/Segment_Delaunay_graph_Linf_traits_2.h \
#    CGAL/Segment_Delaunay_graph_simple_site_2.h \
#    CGAL/Segment_Delaunay_graph_simple_storage_site_2.h \
#    CGAL/Segment_Delaunay_graph_site_2.h \
#    CGAL/Segment_Delaunay_graph_storage_site_2.h \
#    CGAL/Segment_Delaunay_graph_storage_site_with_info_2.h \
#    CGAL/Segment_Delaunay_graph_storage_traits_2.h \
#    CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h \
#    CGAL/Segment_Delaunay_graph_traits_2.h \
#    CGAL/Segment_Delaunay_graph_vertex_base_2.h \
#    CGAL/segment_intersection_points_2.h \
#    CGAL/Segment_tree_d.h \
#    CGAL/Segment_tree_k.h \
#    CGAL/SEP_header.h \
#    CGAL/SEP_to_ImageIO.h \
#    CGAL/Shape_detection_3.h \
#    CGAL/sibson_gradient_fitting.h \
#    CGAL/Side_of_bounded_square_2.h \
#    CGAL/Side_of_oriented_square_2.h \
#    CGAL/Side_of_triangle_mesh.h \
    CGAL/Simple_cartesian.h \
    CGAL/Simple_circular_kernel_2.h \
    CGAL/Simple_homogeneous.h \
    CGAL/Simple_polygon_visibility_2.h \
    CGAL/Simple_spherical_kernel_3.h \
    CGAL/simple_transformations_d.h \
#    CGAL/simplest_rational_in_interval.h \
#    CGAL/Sixtuple.h \
#    CGAL/Skin_surface_3.h \
#    CGAL/Skin_surface_base_3.h \
#    CGAL/Skin_surface_filtered_traits_3.h \
#    CGAL/Skin_surface_marching_tetrahedra_observer_3.h \
#    CGAL/Skin_surface_polyhedral_items_3.h \
#    CGAL/Skin_surface_polyhedral_items_with_face_information.h \
#    CGAL/Skin_surface_quadratic_surface_3.h \
#    CGAL/Skin_surface_refinement_policy_3.h \
#    CGAL/Skin_surface_traits_3.h \
#    CGAL/Skiplist.h \
#    CGAL/Small_side_angle_bisector_decomposition_2.h \
#    CGAL/Snap_rounding_2.h \
#    CGAL/Snap_rounding_kd_2.h \
#    CGAL/Snap_rounding_traits_2.h \
#    CGAL/sorted_matrix_search.h \
#    CGAL/Sorted_matrix_search_traits_adaptor.h \
#    CGAL/Spatial_lock_grid_3.h \
#    CGAL/spatial_sort.h \
#    CGAL/spatial_sort_on_sphere.h \
#    CGAL/Spatial_sort_traits_adapter_2.h \
#    CGAL/Spatial_sort_traits_adapter_3.h \
#    CGAL/Spatial_sort_traits_adapter_d.h \
#    CGAL/Sphere_3.h \
#    CGAL/Spherical_kernel_3.h \
#    CGAL/Spherical_kernel_intersections.h \
#    CGAL/Spherical_kernel_type_equality_wrapper.h \
#    CGAL/Splitters.h \
#    CGAL/Sqrt_extension.h \
#    CGAL/Sqrt_extension_fwd.h \
#    CGAL/squared_distance_2.h \
#    CGAL/squared_distance_2_1.h \
#    CGAL/squared_distance_2_2.h \
#    CGAL/squared_distance_3.h \
#    CGAL/squared_distance_3_0.h \
#    CGAL/squared_distance_3_1.h \
#    CGAL/squared_distance_3_2.h \
#    CGAL/squared_distance_utils.h \
#    CGAL/sse2.h \
#    CGAL/Static_filtered_predicate.h \
#    CGAL/stddef.h \
#    CGAL/Straight_2.h \
#    CGAL/Straight_skeleton_2.h \
#    CGAL/Straight_skeleton_builder_2.h \
#    CGAL/Straight_skeleton_builder_traits_2.h \
#    CGAL/Straight_skeleton_converter_2.h \
#    CGAL/Straight_skeleton_face_base_2.h \
#    CGAL/Straight_skeleton_halfedge_base_2.h \
#    CGAL/Straight_skeleton_items_2.h \
#    CGAL/Straight_skeleton_vertex_base_2.h \
#    CGAL/Stream_lines_2.h \
#    CGAL/streamlines_assertions.h \
#    CGAL/structure_point_set.h \
#    CGAL/subdivide_skin_surface_mesh_3.h \
#    CGAL/subdivide_union_of_balls_mesh_3.h \
#    CGAL/subdivision_method_3.h \
    CGAL/Surface_mesh.h \
    CGAL/Surface_mesh_cell_base_3.h \
    CGAL/Surface_mesh_complex_2_in_triangulation_3.h \
    CGAL/Surface_mesh_default_criteria_3.h \
    CGAL/Surface_mesh_default_edges_criteria_3.h \
    CGAL/Surface_mesh_default_triangulation_3.h \
    CGAL/Surface_mesh_deformation.h \
    CGAL/surface_mesh_parameterization.h \
    CGAL/Surface_mesh_shortest_path.h \
    CGAL/Surface_mesh_traits_generator_3.h \
    CGAL/Surface_mesh_triangulation_generator_3.h \
    CGAL/Surface_mesh_vertex_base_3.h \
#    CGAL/Surface_mesher_generator.h \
#    CGAL/surface_neighbor_coordinates_3.h \
#    CGAL/surface_neighbors_3.h \
#    CGAL/surface_reconstruction_points_assertions.h \
#    CGAL/Surface_sweep_2.h \
#    CGAL/Surface_sweep_2_algorithms.h \
#    CGAL/Sweep_line_2_algorithms.h \
#    CGAL/sweep_observer.h \
#    CGAL/tags.h \
#    CGAL/TDS_full_cell_default_storage_policy.h \
#    CGAL/TDS_full_cell_mirror_storage_policy.h \
#    CGAL/test_FPU_rounding_mode_impl.h \
#    CGAL/Tetrahedron_3.h \
#    CGAL/thread.h \
#    CGAL/Threetuple.h \
#    CGAL/Time_stamper.h \
#    CGAL/Timer.h \
#    CGAL/Timer_impl.h \
#    CGAL/to_rational.h \
#    CGAL/trace.h \
#    CGAL/Transform_iterator.h \
#    CGAL/transforming_iterator.h \
#    CGAL/transforming_pair_iterator.h \
#    CGAL/Tree_assertions.h \
#    CGAL/Tree_base.h \
#    CGAL/Tree_traits.h \
    CGAL/Triangle_2.h \
    CGAL/Triangle_2_Iso_rectangle_2_intersection.h \
    CGAL/Triangle_2_Line_2_intersection.h \
    CGAL/Triangle_2_Point_2_intersection.h \
    CGAL/Triangle_2_Ray_2_intersection.h \
    CGAL/Triangle_2_Segment_2_intersection.h \
    CGAL/Triangle_2_Triangle_2_do_intersect.h \
    CGAL/Triangle_2_Triangle_2_intersection.h \
    CGAL/Triangle_3.h \
    CGAL/Triangle_3_Line_3_do_intersect.h \
    CGAL/Triangle_3_Plane_3_do_intersect.h \
    CGAL/Triangle_3_Point_3_do_intersect.h \
    CGAL/Triangle_3_Ray_3_do_intersect.h \
    CGAL/Triangle_3_Segment_3_do_intersect.h \
    CGAL/Triangle_3_Tetrahedron_3_do_intersect.h \
    CGAL/Triangle_3_Triangle_3_do_intersect.h \
    CGAL/Triangle_3_Triangle_3_intersection.h \
#    CGAL/Triangle_accessor_3.h \
#    CGAL/Triangle_accessor_with_ppmap_3.h \
#    CGAL/Triangular_expansion_visibility_2.h \
#    CGAL/Triangular_field_2.h \
#    CGAL/triangulate_mixed_complex_3.h \
#    CGAL/triangulate_power_diagram_3.h \
#    CGAL/Triangulated_mixed_complex_observer_3.h \
#    CGAL/Triangulation.h \
#    CGAL/Triangulation_2.h \
#    CGAL/Triangulation_2_filtered_projection_traits_3.h \
#    CGAL/Triangulation_2_projection_traits_3.h \
#    CGAL/Triangulation_2_to_lcc.h \
#    CGAL/Triangulation_2_traits_3.h \
#    CGAL/Triangulation_3.h \
#    CGAL/Triangulation_3_to_lcc.h \
#    CGAL/triangulation_assertions.h \
#    CGAL/Triangulation_cell_base_3.h \
#    CGAL/Triangulation_cell_base_with_info_3.h \
#    CGAL/Triangulation_conformer_2.h \
#    CGAL/Triangulation_data_structure.h \
#    CGAL/Triangulation_data_structure_2.h \
#    CGAL/Triangulation_data_structure_3.h \
#    CGAL/Triangulation_data_structure_using_list_2.h \
#    CGAL/Triangulation_default_data_structure_2.h \
#    CGAL/Triangulation_ds_cell_base_3.h \
#    CGAL/Triangulation_ds_circulators_2.h \
#    CGAL/Triangulation_ds_face_2.h \
#    CGAL/Triangulation_ds_face_base_2.h \
#    CGAL/Triangulation_ds_full_cell.h \
#    CGAL/Triangulation_ds_iterators_2.h \
#    CGAL/Triangulation_ds_vertex.h \
#    CGAL/Triangulation_ds_vertex_2.h \
#    CGAL/Triangulation_ds_vertex_base_2.h \
#    CGAL/Triangulation_ds_vertex_base_3.h \
#    CGAL/Triangulation_euclidean_traits_2.h \
#    CGAL/Triangulation_face.h \
#    CGAL/Triangulation_face_base_2.h \
#    CGAL/Triangulation_face_base_with_info_2.h \
#    CGAL/Triangulation_full_cell.h \
#    CGAL/Triangulation_geom_traits_3.h \
#    CGAL/Triangulation_hierarchy_2.h \
#    CGAL/Triangulation_hierarchy_3.h \
#    CGAL/Triangulation_hierarchy_vertex_base_2.h \
#    CGAL/Triangulation_hierarchy_vertex_base_3.h \
#    CGAL/Triangulation_incremental_builder_3.h \
#    CGAL/Triangulation_iterator_adaptator.h \
#    CGAL/Triangulation_line_face_circulator_2.h \
#    CGAL/Triangulation_simplex_3.h \
#    CGAL/Triangulation_structural_filtering_traits.h \
#    CGAL/Triangulation_utils_2.h \
#    CGAL/Triangulation_utils_3.h \
#    CGAL/Triangulation_vertex.h \
#    CGAL/Triangulation_vertex_base_2.h \
#    CGAL/Triangulation_vertex_base_3.h \
#    CGAL/Triangulation_vertex_base_with_id_2.h \
#    CGAL/Triangulation_vertex_base_with_info_2.h \
#    CGAL/Triangulation_vertex_base_with_info_3.h \
#    CGAL/Trivial_iterator.h \
#    CGAL/tss.h \
#    CGAL/tuple.h \
#    CGAL/Twotuple.h \
#    CGAL/type_traits.h \
#    CGAL/typeset.h \
#    CGAL/Umbilics.h \
#    CGAL/Uncertain.h \
#    CGAL/Unfiltered_predicate_adaptor.h \
#    CGAL/Union_find.h \
#    CGAL/Union_of_balls_3.h \
#    CGAL/Unique_hash_map.h \
#    CGAL/unordered.h \
#    CGAL/use.h \
#    CGAL/user_classes.h \
#    CGAL/utility.h \
#    CGAL/utils.h \
#    CGAL/utils_classes.h \
#    CGAL/value_type_traits.h \
#    CGAL/vcm_estimate_edges.h \
#    CGAL/vcm_estimate_normals.h \
    CGAL/vector.h \
    CGAL/Vector_2.h \
    CGAL/Vector_3.h \
#    CGAL/version.h \
#    CGAL/version_macros.h \
#    CGAL/Vertex2Data_Property_Map_with_std_map.h \
#    CGAL/Voronoi_diagram_2.h \
#    CGAL/Voronoi_intersection_2_traits_3.h \
#    CGAL/vtkSurfaceMesherContourFilter.h \
#    CGAL/Weighted_alpha_shape_euclidean_traits_2.h \
#    CGAL/Weighted_alpha_shape_euclidean_traits_3.h \
#    CGAL/Weighted_Minkowski_distance.h \
#    CGAL/Weighted_point.h \
#    CGAL/Weighted_point_2.h \
#    CGAL/Weighted_point_3.h \
#    CGAL/Width_3.h \
#    CGAL/width_assertions.h \
#    CGAL/Width_default_traits_3.h \
#    CGAL/Width_polyhedron_3.h \
#    CGAL/wlop_simplify_and_regularize_point_set.h \
#    CGAL/wmult.h

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../CGAL-4.13.1/build/lib/release/ -lCGAL
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../CGAL-4.13.1/build/lib/debug/ -lCGAL
else:unix: LIBS += -L$$PWD/../../CGAL-4.13.1/build/lib/ -lCGAL

INCLUDEPATH += $$PWD/../../CGAL-4.13.1/build/include
DEPENDPATH += $$PWD/../../CGAL-4.13.1/build/include
