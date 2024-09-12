#pragma once

#include <fstream>

#include <glad/glad.h>
// Glad must be included before GLFW
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <tinyfiledialogs.h>

#include "mesh.h"

class Viewer {
public:
    Viewer(Mesh& mesh);

    void start(std::function<void()> callback, std::string path, std::string title);

    void reset();

    void update();

    bool is_closed();

    void clear_extras();
    void add_triangle(const std::vector<Eigen::Vector3d>& positions, const Eigen::Vector3d& color = {1, 1, 1}, const std::vector<Eigen::Vector3d>& normals = {});
    void add_line(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& color = {1, 1, 1});
    void add_sphere(const Eigen::Vector3d& c, double r, const Eigen::Vector3d& color = {1, 1, 1}, int res = 2);
    void add_cylinder(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double r, const Eigen::Vector3d& color = {1, 1, 1}, int res = 8);
    void add_cone(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double r, const Eigen::Vector3d& color = {1, 1, 1}, int res = 16);
    void add_coordinate_system();

    void set_scale(double scale);

    bool show_cells_;
    bool show_faces_;
    bool show_edges_;
    bool show_vertices_;

    bool vertex_face_colors_;
    bool use_texture_;
    bool parameter_space_;
    bool twod_;
    std::function<void(int type, int idx, const Eigen::Vector3d& p)> pick_callback_;
    bool show_boundary_;
    bool show_colored_;
    bool rotate_;
    bool paper_;

private:
    void frame();

    Eigen::Vector3d position(Vertex v);

    void to_float_mat(const Eigen::Matrix4d& m, float array[16]);

    Eigen::Vector2d to_screen(ImVec2 mouse_pos, ImVec2 size);

    Mesh& mesh_;

    std::function<void()> callback_;

    GLFWwindow* window_;
    double pixel_ratio_;
    unsigned int shader_program_;
    unsigned int texture_;
    unsigned int framebuffer_;
    unsigned int colorbuffer_;
    unsigned int depthbuffer_;
    unsigned int cell_vao_;
    unsigned int cell_position_buffer_;
    unsigned int cell_normal_buffer_;
    unsigned int cell_color_buffer_;
    unsigned int cell_centroid_buffer_;
    unsigned int cell_index_buffer_;
    unsigned int n_cell_indices_;
    unsigned int face_vao_;
    unsigned int face_position_buffer_;
    unsigned int face_normal_buffer_;
    unsigned int face_color_buffer_;
    unsigned int tex_coord_buffer_;
    unsigned int face_centroid_buffer_;
    unsigned int face_index_buffer_;
    unsigned int n_face_indices_;
    unsigned int edge_vao_;
    unsigned int edge_position_buffer_;
    unsigned int edge_color_buffer_;
    unsigned int edge_index_buffer_;
    unsigned int n_edge_indices_;
    unsigned int vertex_instance_vao_;
    unsigned int vertex_instance_position_buffer_;
    unsigned int vertex_instance_normal_buffer_;
    unsigned int vertex_instance_index_buffer_;
    unsigned int n_vertex_instance_indices_;
    unsigned int vertex_position_buffer_;
    unsigned int vertex_color_buffer_;
    unsigned int n_vertex_instances_;
    unsigned int extra_vao_;
    unsigned int extra_position_buffer_;
    unsigned int extra_normal_buffer_;
    unsigned int extra_color_buffer_;
    unsigned int extra_index_buffer_;
    unsigned int n_extra_indices_;
    unsigned int extra_line_index_buffer_;
    unsigned int n_extra_line_indices_;

    std::chrono::steady_clock::time_point prev_frame_;
    double size_;
    Eigen::Matrix4d transformation_matrix_;
    Eigen::Vector2d m0_;
    bool rot_;
    bool trans_;
    bool rot_cam_;
    double scale_;
    bool closed_;

    std::vector<std::vector<Eigen::Vector3d>> extra_triangle_positions_;
    std::vector<std::vector<Eigen::Vector3d>> extra_triangle_normals_;
    std::vector<std::vector<Eigen::Vector3d>> extra_triangle_colors_;
    std::vector<std::vector<Eigen::Vector3d>> extra_line_positions_;
    std::vector<std::vector<Eigen::Vector3d>> extra_line_colors_;
};
