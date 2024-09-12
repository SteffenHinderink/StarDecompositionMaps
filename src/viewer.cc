#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#include "viewer.h"

Viewer::Viewer(Mesh& mesh)
    : mesh_(mesh), rot_(false), trans_(false), rot_cam_(false), scale_(0.5), closed_(false),
      show_cells_(true), show_faces_(false), show_edges_(false), show_vertices_(false),
      vertex_face_colors_(false), use_texture_(false), parameter_space_(false), twod_(false),
      pick_callback_([](int type, int idx, const Eigen::Vector3d& p) {}),
      show_boundary_(false), show_colored_(false), rotate_(false), paper_(false) {
#ifndef TET
    scale_ = 1;
    show_faces_ = true;
    show_cells_ = false;
#endif
}

void Viewer::start(std::function<void()> callback, std::string path, std::string title) {
    callback_ = callback;
    glfwSetErrorCallback([](int error, const char* description) {
        std::cout << "GLFW Error " << error << ": " << description << std::endl;
    });
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    window_ = glfwCreateWindow(1920, 1080, title.c_str(), nullptr, nullptr);
    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1);
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::GetIO().ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGui_ImplGlfw_InitForOpenGL(window_, true);
    ImGui_ImplOpenGL3_Init("#version 410");
    gladLoadGLLoader((GLADloadproc) glfwGetProcAddress);

    int wx, fx, tmp;
    glfwGetWindowSize(window_, &wx, &tmp);
    glfwGetFramebufferSize(window_, &fx, &tmp);
    pixel_ratio_ = (double) fx / wx;

    int success;
    char info_log[1024];
    std::ifstream file;

    file = std::ifstream(path + SEP + "shader" + SEP + "vs.vert");
    std::string vs((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    const char* vertex_shader_source = vs.c_str();
    unsigned int vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, &vertex_shader_source, nullptr);
    glCompileShader(vertex_shader);
    glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertex_shader, 1024, NULL, info_log);
        std::cout << "Vertex shader compilation failed" << std::endl;
        std::cout << info_log << std::endl;
    }

    file = std::ifstream(path + SEP + "shader" + SEP + "fs.frag");
    std::string fs((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    const char* fragment_shader_source = fs.c_str();
    unsigned int fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &fragment_shader_source, nullptr);
    glCompileShader(fragment_shader);
    glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragment_shader, 1024, NULL, info_log);
        std::cout << "Fragment shader compilation failed" << std::endl;
        std::cout << info_log << std::endl;
    }

    shader_program_ = glCreateProgram();
    glAttachShader(shader_program_, vertex_shader);
    glAttachShader(shader_program_, fragment_shader);
    glLinkProgram(shader_program_);
    glUseProgram(shader_program_);
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);

    glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    stbi_set_flip_vertically_on_load(true);
    int tx, ty;
    unsigned char* texture_data = stbi_load((path + SEP + "res" + SEP + "texture.png").c_str(), &tx, &ty, &tmp, 0);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tx, ty, 0, GL_RGB, GL_UNSIGNED_BYTE, texture_data);
    glGenerateMipmap(GL_TEXTURE_2D);
    stbi_image_free(texture_data);

    glGenFramebuffers(1, &framebuffer_);
    glGenTextures(1, &colorbuffer_);
    glGenRenderbuffers(1, &depthbuffer_);
    glGenVertexArrays(1, &cell_vao_);
    glGenBuffers(1, &cell_position_buffer_);
    glGenBuffers(1, &cell_normal_buffer_);
    glGenBuffers(1, &cell_color_buffer_);
    glGenBuffers(1, &cell_centroid_buffer_);
    glGenBuffers(1, &cell_index_buffer_);
    glGenVertexArrays(1, &face_vao_);
    glGenBuffers(1, &face_position_buffer_);
    glGenBuffers(1, &face_normal_buffer_);
    glGenBuffers(1, &face_color_buffer_);
    glGenBuffers(1, &tex_coord_buffer_);
    glGenBuffers(1, &face_centroid_buffer_);
    glGenBuffers(1, &face_index_buffer_);
    glGenVertexArrays(1, &edge_vao_);
    glGenBuffers(1, &edge_position_buffer_);
    glGenBuffers(1, &edge_color_buffer_);
    glGenBuffers(1, &edge_index_buffer_);
    glGenVertexArrays(1, &vertex_instance_vao_);
    glGenBuffers(1, &vertex_instance_position_buffer_);
    glGenBuffers(1, &vertex_instance_normal_buffer_);
    glGenBuffers(1, &vertex_instance_index_buffer_);
    glGenBuffers(1, &vertex_position_buffer_);
    glGenBuffers(1, &vertex_color_buffer_);
    glGenVertexArrays(1, &extra_vao_);
    glGenBuffers(1, &extra_position_buffer_);
    glGenBuffers(1, &extra_normal_buffer_);
    glGenBuffers(1, &extra_color_buffer_);
    glGenBuffers(1, &extra_index_buffer_);
    glGenBuffers(1, &extra_line_index_buffer_);

    reset();

#ifndef _WIN32
    glfwSetWindowUserPointer(window_, this);
    glfwSetWindowRefreshCallback(window_, [](GLFWwindow* window) {
        ((Viewer*) glfwGetWindowUserPointer(window))->frame();
    });
#endif

    while (!glfwWindowShouldClose(window_)) {
        frame();
    }

    glDeleteProgram(shader_program_);
    glDeleteTextures(1, &texture_);
    glDeleteFramebuffers(1, &framebuffer_);
    glDeleteTextures(1, &colorbuffer_);
    glDeleteRenderbuffers(1, &depthbuffer_);
    glDeleteVertexArrays(1, &cell_vao_);
    glDeleteBuffers(1, &cell_position_buffer_);
    glDeleteBuffers(1, &cell_normal_buffer_);
    glDeleteBuffers(1, &cell_color_buffer_);
    glDeleteBuffers(1, &cell_centroid_buffer_);
    glDeleteBuffers(1, &cell_index_buffer_);
    glDeleteVertexArrays(1, &face_vao_);
    glDeleteBuffers(1, &face_position_buffer_);
    glDeleteBuffers(1, &face_normal_buffer_);
    glDeleteBuffers(1, &face_color_buffer_);
    glDeleteBuffers(1, &tex_coord_buffer_);
    glDeleteBuffers(1, &face_centroid_buffer_);
    glDeleteBuffers(1, &face_index_buffer_);
    glDeleteVertexArrays(1, &edge_vao_);
    glDeleteBuffers(1, &edge_position_buffer_);
    glDeleteBuffers(1, &edge_color_buffer_);
    glDeleteBuffers(1, &edge_index_buffer_);
    glDeleteVertexArrays(1, &vertex_instance_vao_);
    glDeleteBuffers(1, &vertex_instance_position_buffer_);
    glDeleteBuffers(1, &vertex_instance_normal_buffer_);
    glDeleteBuffers(1, &vertex_instance_index_buffer_);
    glDeleteBuffers(1, &vertex_position_buffer_);
    glDeleteBuffers(1, &vertex_color_buffer_);
    glDeleteVertexArrays(1, &extra_vao_);
    glDeleteBuffers(1, &extra_position_buffer_);
    glDeleteBuffers(1, &extra_normal_buffer_);
    glDeleteBuffers(1, &extra_color_buffer_);
    glDeleteBuffers(1, &extra_index_buffer_);
    glDeleteBuffers(1, &extra_line_index_buffer_);

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window_);
    glfwTerminate();

    closed_ = true;
}

void Viewer::reset() {
    Eigen::Vector3d min = Eigen::Vector3d::Constant(std::numeric_limits<double>::max());
    Eigen::Vector3d max = Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest());
    for (auto v : mesh_.vertices()) {
        Eigen::Vector3d p = position(v);
        for (int i = 0; i < 3; i++) {
            if (p(i) < min(i)) {
                min(i) = p(i);
            }
            if (p(i) > max(i)) {
                max(i) = p(i);
            }
        }
    }
    size_ = (max - min).norm();
    if (size_ == 0) {
        size_ = 1;
    }
    Eigen::Vector3d c = (min + max) / 2;
    transformation_matrix_ = Eigen::Matrix4d::Identity();
    transformation_matrix_.block<3, 1>(0, 3) = -c;

    update();
}

void Viewer::update() {
    int position_location = glGetAttribLocation(shader_program_, "aPosition");
    int normal_location = glGetAttribLocation(shader_program_, "aNormal");
    int color_location = glGetAttribLocation(shader_program_, "aColor");
    int tex_coord_location = glGetAttribLocation(shader_program_, "aTexCoord");
    int centroid_location = glGetAttribLocation(shader_program_, "aCentroid");
    int instanced_position_location = glGetAttribLocation(shader_program_, "aInstancedPosition");

    Eigen::Vector3d sel_color = {0, 1, 0.8};
#ifdef TET
    {
        std::vector<float> positions;
        std::vector<float> normals;
        std::vector<float> colors;
        std::vector<float> centroids;
        std::vector<int> indices;

        auto sel = mesh_.property<Cell, bool>("sel", false);
        auto clr = mesh_.property<Cell, Eigen::Vector3d>("clr", {1, 1, 1});
        for (auto c : mesh_.cells()) {
            if (show_boundary_ && !mesh_.is_boundary(c)) {
                continue;
            }
            Eigen::Vector3d color = clr[c];
            if ((show_colored_ && color == Eigen::Vector3d(1, 1, 1)) || color == Eigen::Vector3d(-1, -1, -1)) {
                continue;
            }
            if (sel[c]) {
                color = sel_color;
            }
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
            for (auto cv : mesh_.tet_vertices(c)) {
                centroid += position(cv);
            }
            centroid /= 4;
            for (auto ch : mesh_.cell_halffaces(c)) {
                std::vector<Eigen::Vector3d> halfface_positions;
                for (auto chv : mesh_.halfface_vertices(ch)) {
                    halfface_positions.push_back(position(chv));
                }
                std::vector<Eigen::Vector3d> halfface_normals = std::vector<Eigen::Vector3d>(3, (halfface_positions[1] - halfface_positions[0]).cross(halfface_positions[2] - halfface_positions[0]).normalized());
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        positions.push_back(halfface_positions[i](j));
                        normals.push_back(halfface_normals[i](j));
                        colors.push_back(color(j));
                        centroids.push_back(centroid(j));
                    }
                    indices.push_back(indices.size());
                }
            }
        }
        n_cell_indices_ = indices.size();

        glBindVertexArray(cell_vao_);

        glBindBuffer(GL_ARRAY_BUFFER, cell_position_buffer_);
        glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(float), positions.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(position_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(position_location);

        glBindBuffer(GL_ARRAY_BUFFER, cell_normal_buffer_);
        glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), normals.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(normal_location, 3, GL_FLOAT, GL_TRUE, 0, 0);
        glEnableVertexAttribArray(normal_location);

        glBindBuffer(GL_ARRAY_BUFFER, cell_color_buffer_);
        glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(float), colors.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(color_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(color_location);

        glBindBuffer(GL_ARRAY_BUFFER, cell_centroid_buffer_);
        glBufferData(GL_ARRAY_BUFFER, centroids.size() * sizeof(float), centroids.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(centroid_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(centroid_location);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell_index_buffer_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_cell_indices_ * sizeof(int), indices.data(), GL_STATIC_DRAW);
    }
#endif
    {
        std::vector<float> positions;
        std::vector<float> normals;
        std::vector<float> colors;
        std::vector<float> tex_coords;
        std::vector<float> centroids;
        std::vector<int> indices;

        auto sel = mesh_.property<Face, bool>("sel", false);
        auto clr = mesh_.property<Face, Eigen::Vector3d>("clr", {1, 1, 1});
        auto vertex_clr = mesh_.property<Vertex, Eigen::Vector3d>("clr", {1, 1, 1});
        auto tex = mesh_.property<Vertex, Eigen::Vector2d>("tex", {0.5, 0.5});
        for (auto f : mesh_.faces()) {
            if (show_boundary_ && !mesh_.is_boundary(f)) {
                continue;
            }
            Eigen::Vector3d color = clr[f];
            if ((show_colored_ && color == Eigen::Vector3d(1, 1, 1)) || color == Eigen::Vector3d(-1, -1, -1)) {
                continue;
            }
            if (sel[f]) {
                color = sel_color;
            }
            std::vector<Eigen::Vector3d> face_positions;
            std::vector<Eigen::Vector3d> face_colors = std::vector<Eigen::Vector3d>(3, color);
            std::vector<Eigen::Vector2d> face_tex_coords;
            Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
            int i = 0;
            for (auto fv : mesh_.face_vertices(f)) {
                face_positions.push_back(position(fv));
                if (vertex_face_colors_ && !sel[f]) {
                    face_colors[i] = vertex_clr[fv];
                }
                face_tex_coords.push_back(tex[fv]);
                centroid += position(fv);
                i++;
            }
            centroid /= 3;
            std::vector<Eigen::Vector3d> face_normals = std::vector<Eigen::Vector3d>(3, (face_positions[1] - face_positions[0]).cross(face_positions[2] - face_positions[0]).normalized());
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    positions.push_back(face_positions[i](j));
                    normals.push_back(face_normals[i](j));
                    colors.push_back(face_colors[i](j));
                    centroids.push_back(centroid(j));
                }
                for (int j = 0; j < 2; j++) {
                    tex_coords.push_back(face_tex_coords[i](j));
                }
                indices.push_back(indices.size());
            }
        }
        n_face_indices_ = indices.size();

        glBindVertexArray(face_vao_);

        glBindBuffer(GL_ARRAY_BUFFER, face_position_buffer_);
        glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(float), positions.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(position_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(position_location);

        glBindBuffer(GL_ARRAY_BUFFER, face_normal_buffer_);
        glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), normals.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(normal_location, 3, GL_FLOAT, GL_TRUE, 0, 0);
        glEnableVertexAttribArray(normal_location);

        glBindBuffer(GL_ARRAY_BUFFER, face_color_buffer_);
        glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(float), colors.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(color_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(color_location);

        glBindBuffer(GL_ARRAY_BUFFER, tex_coord_buffer_);
        glBufferData(GL_ARRAY_BUFFER, tex_coords.size() * sizeof(float), tex_coords.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(tex_coord_location, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(tex_coord_location);

        glBindBuffer(GL_ARRAY_BUFFER, face_centroid_buffer_);
        glBufferData(GL_ARRAY_BUFFER, centroids.size() * sizeof(float), centroids.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(centroid_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(centroid_location);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, face_index_buffer_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_face_indices_ * sizeof(int), indices.data(), GL_STATIC_DRAW);
    }
    {
        std::vector<float> positions;
        std::vector<float> colors;
        std::vector<int> indices;

        auto sel = mesh_.property<Edge, bool>("sel", false);
        auto clr = mesh_.property<Edge, Eigen::Vector3d>("clr", {1, 1, 1});
        for (auto e : mesh_.edges()) {
            if (show_boundary_ && !mesh_.is_boundary(e)) {
                continue;
            }
            Eigen::Vector3d color = clr[e];
            if ((show_colored_ && color == Eigen::Vector3d(1, 1, 1)) || color == Eigen::Vector3d(-1, -1, -1)) {
                continue;
            }
            if (sel[e]) {
                color = sel_color;
            }
            Halfedge h = mesh_.halfedge_handle(e, 0);
            std::vector<Eigen::Vector3d> edge_positions = {
                position(mesh_.from_vertex_handle(h)),
                position(mesh_.to_vertex_handle(h))};
            std::vector<Eigen::Vector3d> edge_colors = std::vector<Eigen::Vector3d>(2, color);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 3; j++) {
                    positions.push_back(edge_positions[i](j));
                    colors.push_back(edge_colors[i](j));
                }
                indices.push_back(indices.size());
            }
        }
        n_edge_indices_ = indices.size();

        glBindVertexArray(edge_vao_);

        glBindBuffer(GL_ARRAY_BUFFER, edge_position_buffer_);
        glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(float), positions.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(position_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(position_location);
        glVertexAttribPointer(centroid_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(centroid_location);

        glBindBuffer(GL_ARRAY_BUFFER, edge_color_buffer_);
        glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(float), colors.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(color_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(color_location);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_index_buffer_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_edge_indices_ * sizeof(int), indices.data(), GL_STATIC_DRAW);
    }
    {
        std::vector<float> instance_positions;
        std::vector<float> instance_normals;
        std::vector<int> instance_indices;
        std::vector<float> positions;
        std::vector<float> colors;

        auto add_triangle = [&](const std::vector<Eigen::Vector3d>& positions, const std::vector<Eigen::Vector3d>& normals) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    instance_positions.push_back(positions[i](j));
                    instance_normals.push_back(normals[i](j));
                }
                instance_indices.push_back(instance_indices.size());
            }
        };
        Eigen::Vector3d c = Eigen::Vector3d::Zero();
        double r = size_ / 1024;
        int res = 2;
        double phi = 1.61803398875; // = (1 + sqrt(5)) / 2 (Golden ratio)
        std::vector<Eigen::Vector3d> vertices = {
            {0, 1, phi},
            {0, 1, -phi},
            {0, -1, phi},
            {0, -1, -phi},
            {1, phi, 0},
            {1, -phi, 0},
            {-1, phi, 0},
            {-1, -phi, 0},
            {phi, 0, 1},
            {-phi, 0, 1},
            {phi, 0, -1},
            {-phi, 0, -1}};
        std::vector<std::vector<int>> faces = {
            {0, 4, 6},
            {4, 1, 6},
            {6, 1, 11},
            {1, 3, 11},
            {11, 3, 7},
            {3, 5, 7},
            {7, 5, 2},
            {5, 8, 2},
            {2, 8, 0},
            {8, 4, 0},
            {9, 0, 6},
            {9, 6, 11},
            {9, 11, 7},
            {9, 7, 2},
            {9, 2, 0},
            {10, 4, 8},
            {10, 8, 5},
            {10, 5, 3},
            {10, 3, 1},
            {10, 1, 4}};
        for (int i = 0; i < faces.size(); i++) {
            std::vector<Eigen::Vector3d> triangle(3);
            for (int j = 0; j < 3; j++) {
                triangle[j] = vertices[faces[i][j]];
            }
            Eigen::Vector3d ab = (triangle[1] - triangle[0]) / res;
            Eigen::Vector3d ac = (triangle[2] - triangle[0]) / res;
            std::vector<std::vector<Eigen::Vector3d>> p(res + 1);
            std::vector<std::vector<Eigen::Vector3d>> n(res + 1);
            for (int j = 0; j <= res; j++) {
                p[j] = std::vector<Eigen::Vector3d>(res - j + 1);
                n[j] = std::vector<Eigen::Vector3d>(res - j + 1);
                for (int k = 0; k <= res - j; k++) {
                    n[j][k] = (triangle[0] + j * ab + k * ac).normalized();
                    p[j][k] = c + r * n[j][k];
                }
            }
            for (int j = 0; j < res; j++) {
                for (int k = 0; k < res - j; k++) {
                    add_triangle({p[j][k], p[j + 1][k], p[j][k + 1]}, {n[j][k], n[j + 1][k], n[j][k + 1]});
                    if (k < res - j - 1) {
                        add_triangle({p[j][k + 1], p[j + 1][k], p[j + 1][k + 1]}, {n[j][k + 1], n[j + 1][k], n[j + 1][k + 1]});
                    }
                }
            }
        }
        n_vertex_instance_indices_ = instance_indices.size();
        auto sel = mesh_.property<Vertex, bool>("sel", false);
        auto clr = mesh_.property<Vertex, Eigen::Vector3d>("clr", {1, 1, 1});
        for (auto v : mesh_.vertices()) {
            if (show_boundary_ && !mesh_.is_boundary(v)) {
                continue;
            }
            Eigen::Vector3d color = clr[v];
            if ((show_colored_ && color == Eigen::Vector3d(1, 1, 1)) || color == Eigen::Vector3d(-1, -1, -1)) {
                continue;
            }
            if (sel[v]) {
                color = sel_color;
            }
            Eigen::Vector3d vertex_position = position(v);
            for (int i = 0; i < 3; i++) {
                positions.push_back(vertex_position(i));
                colors.push_back(color(i));
            }
        }
        n_vertex_instances_ = positions.size() / 3;

        glBindVertexArray(vertex_instance_vao_);

        glBindBuffer(GL_ARRAY_BUFFER, vertex_instance_position_buffer_);
        glBufferData(GL_ARRAY_BUFFER, instance_positions.size() * sizeof(float), instance_positions.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(position_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(position_location);

        glBindBuffer(GL_ARRAY_BUFFER, vertex_instance_normal_buffer_);
        glBufferData(GL_ARRAY_BUFFER, instance_normals.size() * sizeof(float), instance_normals.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(normal_location, 3, GL_FLOAT, GL_TRUE, 0, 0);
        glEnableVertexAttribArray(normal_location);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vertex_instance_index_buffer_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_vertex_instance_indices_ * sizeof(int), instance_indices.data(), GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, vertex_position_buffer_);
        glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(float), positions.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(instanced_position_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glVertexAttribPointer(centroid_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glVertexAttribDivisor(instanced_position_location, 1);
        glEnableVertexAttribArray(instanced_position_location);
        glVertexAttribDivisor(centroid_location, 1);
        glEnableVertexAttribArray(centroid_location);

        glBindBuffer(GL_ARRAY_BUFFER, vertex_color_buffer_);
        glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(float), colors.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(color_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glVertexAttribDivisor(color_location, 1);
        glEnableVertexAttribArray(color_location);
    }
    {
        std::vector<float> extra_positions;
        std::vector<float> extra_normals;
        std::vector<float> extra_colors;
        std::vector<int> extra_indices;
        std::vector<int> extra_line_indices;

        for (int i = 0; i < extra_triangle_positions_.size(); i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    extra_positions.push_back(extra_triangle_positions_[i][j](k));
                    extra_normals.push_back(extra_triangle_normals_[i][j](k));
                    extra_colors.push_back(extra_triangle_colors_[i][j](k));
                }
                extra_indices.push_back(extra_indices.size());
            }
        }
        n_extra_indices_ = extra_indices.size();
        for (int i = 0; i < extra_line_positions_.size(); i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 3; k++) {
                    extra_positions.push_back(extra_line_positions_[i][j](k));
                    extra_colors.push_back(extra_line_colors_[i][j](k));
                }
                extra_line_indices.push_back(n_extra_indices_ + extra_line_indices.size());
            }
        }
        n_extra_line_indices_ = extra_line_indices.size();

        glBindVertexArray(extra_vao_);

        glBindBuffer(GL_ARRAY_BUFFER, extra_position_buffer_);
        glBufferData(GL_ARRAY_BUFFER, extra_positions.size() * sizeof(float), extra_positions.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(position_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(position_location);
        glVertexAttribPointer(centroid_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(centroid_location);

        if (n_extra_indices_ > 0) {
            glBindBuffer(GL_ARRAY_BUFFER, extra_normal_buffer_);
            glBufferData(GL_ARRAY_BUFFER, extra_normals.size() * sizeof(float), extra_normals.data(), GL_STATIC_DRAW);
            glVertexAttribPointer(normal_location, 3, GL_FLOAT, GL_TRUE, 0, 0);
            glEnableVertexAttribArray(normal_location);
        }

        glBindBuffer(GL_ARRAY_BUFFER, extra_color_buffer_);
        glBufferData(GL_ARRAY_BUFFER, extra_colors.size() * sizeof(float), extra_colors.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(color_location, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(color_location);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, extra_index_buffer_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_extra_indices_ * sizeof(int), extra_indices.data(), GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, extra_line_index_buffer_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_extra_line_indices_ * sizeof(int), extra_line_indices.data(), GL_STATIC_DRAW);
    }
}

bool Viewer::is_closed() {
    return closed_;
}

void Viewer::clear_extras() {
    extra_triangle_positions_.clear();
    extra_triangle_normals_.clear();
    extra_triangle_colors_.clear();
    extra_line_positions_.clear();
    extra_line_colors_.clear();
}

void Viewer::add_triangle(const std::vector<Eigen::Vector3d>& positions, const Eigen::Vector3d& color, const std::vector<Eigen::Vector3d>& _normals) {
    std::vector<Eigen::Vector3d> normals = _normals;
    if (normals.size() == 0) {
        normals = std::vector<Eigen::Vector3d>(3, (positions[1] - positions[0]).cross(positions[2] - positions[0]).normalized());
    }
    extra_triangle_positions_.push_back(positions);
    extra_triangle_normals_.push_back(normals);
    extra_triangle_colors_.push_back(std::vector<Eigen::Vector3d>(3, color));
}

void Viewer::add_line(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& color) {
    extra_line_positions_.push_back({a, b});
    extra_line_colors_.push_back(std::vector<Eigen::Vector3d>(2, color));
}

void Viewer::add_sphere(const Eigen::Vector3d& c, double r, const Eigen::Vector3d& color, int res) {
    double phi = 1.61803398875; // = (1 + sqrt(5)) / 2 (Golden ratio)
    std::vector<Eigen::Vector3d> vertices = {
        {0, 1, phi},
        {0, 1, -phi},
        {0, -1, phi},
        {0, -1, -phi},
        {1, phi, 0},
        {1, -phi, 0},
        {-1, phi, 0},
        {-1, -phi, 0},
        {phi, 0, 1},
        {-phi, 0, 1},
        {phi, 0, -1},
        {-phi, 0, -1}};
    std::vector<std::vector<int>> faces = {
        {0, 4, 6},
        {4, 1, 6},
        {6, 1, 11},
        {1, 3, 11},
        {11, 3, 7},
        {3, 5, 7},
        {7, 5, 2},
        {5, 8, 2},
        {2, 8, 0},
        {8, 4, 0},
        {9, 0, 6},
        {9, 6, 11},
        {9, 11, 7},
        {9, 7, 2},
        {9, 2, 0},
        {10, 4, 8},
        {10, 8, 5},
        {10, 5, 3},
        {10, 3, 1},
        {10, 1, 4}};
    for (int i = 0; i < faces.size(); i++) {
        std::vector<Eigen::Vector3d> triangle(3);
        for (int j = 0; j < 3; j++) {
            triangle[j] = vertices[faces[i][j]];
        }
        Eigen::Vector3d ab = (triangle[1] - triangle[0]) / res;
        Eigen::Vector3d ac = (triangle[2] - triangle[0]) / res;
        std::vector<std::vector<Eigen::Vector3d>> p(res + 1);
        std::vector<std::vector<Eigen::Vector3d>> n(res + 1);
        for (int j = 0; j <= res; j++) {
            p[j] = std::vector<Eigen::Vector3d>(res - j + 1);
            n[j] = std::vector<Eigen::Vector3d>(res - j + 1);
            for (int k = 0; k <= res - j; k++) {
                n[j][k] = (triangle[0] + j * ab + k * ac).normalized();
                p[j][k] = c + r * n[j][k];
            }
        }
        for (int j = 0; j < res; j++) {
            for (int k = 0; k < res - j; k++) {
                add_triangle({p[j][k], p[j + 1][k], p[j][k + 1]}, color, {n[j][k], n[j + 1][k], n[j][k + 1]});
                if (k < res - j - 1) {
                    add_triangle({p[j][k + 1], p[j + 1][k], p[j + 1][k + 1]}, color, {n[j][k + 1], n[j + 1][k], n[j + 1][k + 1]});
                }
            }
        }
    }
}

void Viewer::add_cylinder(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double r, const Eigen::Vector3d& color, int res) {
    Eigen::Vector3d x = (b - a).normalized();
    Eigen::Vector3d y = x.cross(std::abs(x(0)) > 0.999 ? Eigen::Vector3d(0, 1, 0) : Eigen::Vector3d(1, 0, 0)).normalized();
    Eigen::Vector3d z = x.cross(y).normalized();
    Eigen::Vector3d n0 = z;
    Eigen::Vector3d a0 = a + r * n0;
    Eigen::Vector3d b0 = b + r * n0;
    for (int i = 0; i < res; i++) {
        double angle = 2 * M_PI * (i + 1) / res;
        Eigen::Vector3d n1 = std::sin(angle) * y + std::cos(angle) * z;
        Eigen::Vector3d a1 = a + r * n1;
        Eigen::Vector3d b1 = b + r * n1;
        add_triangle({a, a0, a1}, color);
        add_triangle({a0, b0, b1}, color, {n0, n0, n1});
        add_triangle({a0, b1, a1}, color, {n0, n1, n1});
        add_triangle({b, b1, b0}, color);
        n0 = n1;
        a0 = a1;
        b0 = b1;
    }
}

void Viewer::add_cone(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double r, const Eigen::Vector3d& color, int res) {
    Eigen::Vector3d ab = b - a;
    double l = ab.norm();
    Eigen::Vector3d x = ab / l;
    Eigen::Vector3d y = x.cross(std::abs(x(0)) > 0.999 ? Eigen::Vector3d(0, 1, 0) : Eigen::Vector3d(1, 0, 0)).normalized();
    Eigen::Vector3d z = x.cross(y).normalized();
    Eigen::Vector3d n0 = (l * z + r * x).normalized();
    Eigen::Vector3d a0 = a + r * z;
    Eigen::Vector3d b0 = b + r * z;
    for (int i = 0; i < res; i++) {
        double half_angle = 2 * M_PI * (i + 0.5) / res;
        Eigen::Vector3d n = (l * (std::sin(half_angle) * y + std::cos(half_angle) * z) + r * x).normalized();
        double angle = 2 * M_PI * (i + 1) / res;
        Eigen::Vector3d d = std::sin(angle) * y + std::cos(angle) * z;
        Eigen::Vector3d n1 = (l * d + r * x).normalized();
        Eigen::Vector3d a1 = a + r * d;
        Eigen::Vector3d b1 = b + r * d;
        add_triangle({a, a0, a1}, color);
        add_triangle({a0, b, a1}, color, {n0, n, n1});
        n0 = n1;
        a0 = a1;
        b0 = b1;
    }
}

void Viewer::add_coordinate_system() {
    for (int i = 0; i < 3; i++) {
        Eigen::Vector3d unit = Eigen::Vector3d::Zero();
        unit(i) = 1;
        add_line({0, 0, 0}, unit, unit);
    }
}

void Viewer::set_scale(double scale) {
    scale_ = scale;
}

void Viewer::frame() {
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    std::chrono::steady_clock::time_point curr_frame = std::chrono::steady_clock::now();
    double time = (double) std::chrono::duration_cast<std::chrono::microseconds>(curr_frame - prev_frame_).count() / 1000000;

    ImGui::Begin("Menu");
    if (ImGui::TreeNode("Mesh info")) {
#ifdef TET
        ImGui::Text("Number of vertices: %lu", mesh_.n_logical_vertices());
        ImGui::Text("Number of halfedges: %lu", mesh_.n_logical_halfedges());
        ImGui::Text("Number of edges: %lu", mesh_.n_logical_edges());
        ImGui::Text("Number of halffaces: %lu", mesh_.n_logical_halffaces());
        ImGui::Text("Number of faces: %lu", mesh_.n_logical_faces());
        ImGui::Text("Number of cells: %lu", mesh_.n_logical_cells());
#else
        ImGui::Text("Number of vertices: %lu", mesh_.n_vertices());
        ImGui::Text("Number of halfedges: %lu", mesh_.n_halfedges());
        ImGui::Text("Number of edges: %lu", mesh_.n_edges());
        ImGui::Text("Number of faces: %lu", mesh_.n_faces());
#endif
        ImGui::TreePop();
        ImGui::Spacing();
    }
#ifdef TET
    ImGui::Checkbox("Show cells", &show_cells_);
    if (ImGui::IsKeyPressed(ImGuiKey_C)) {
        show_cells_ = !show_cells_;
    }
#endif
    ImGui::Checkbox("Show faces", &show_faces_);
    if (ImGui::IsKeyPressed(ImGuiKey_F)) {
        show_faces_ = !show_faces_;
    }
    ImGui::Checkbox("Show edges", &show_edges_);
    if (ImGui::IsKeyPressed(ImGuiKey_E)) {
        show_edges_ = !show_edges_;
    }
    ImGui::Checkbox("Show vertices", &show_vertices_);
    if (ImGui::IsKeyPressed(ImGuiKey_V)) {
        show_vertices_ = !show_vertices_;
    }
    callback_();
    ImGui::End();

    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowSize(io.DisplaySize);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("_", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoBringToFrontOnFocus);
    ImVec2 size = ImGui::GetContentRegionAvail();

    // Handle input
    double fov = 70;
    double fov_factor = 1 / std::tan(fov / 2 * M_PI / 180);
    Eigen::Matrix4d look_at_matrix = Eigen::Matrix4d::Identity();
    look_at_matrix(2, 3) = -size_;
    if (ImGui::IsWindowFocused() && ImGui::IsMouseClicked(ImGuiMouseButton_Left) && (ImGui::IsKeyDown(ImGuiKey_LeftAlt) || twod_)) {
        // Ray
        Eigen::Matrix4d model_view_matrix_inverse = (look_at_matrix * transformation_matrix_).inverse();
        Eigen::Vector2d m = to_screen(ImGui::GetMousePos(), size);
        Eigen::Vector3d target_cam;
        target_cam << m, -(size.x > size.y ? size.y : size.x) * fov_factor / 2;
        Eigen::Vector3d target = (model_view_matrix_inverse * target_cam.homogeneous()).hnormalized();
        Eigen::Vector3d o = (model_view_matrix_inverse * Eigen::Vector4d{0, 0, 0, 1}).hnormalized();
        Eigen::Vector3d d = (target - o).normalized();
        if (twod_) {
            double alpha = -target(2) / (o(2) - target(2));
            o = alpha * o + (1 - alpha) * target;
            o(2) = size_;
            d = {0, 0, -1};
        }
        // Cast
        double min = std::numeric_limits<double>::max();
        int type = -1;
        int idx = -1;
#ifdef TET
        if (show_cells_ || twod_) {
            for (auto c : mesh_.cells()) {
                Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
                for (auto cv : mesh_.tet_vertices(c)) {
                    centroid += position(cv);
                }
                centroid /= 4;
                for (auto ch : mesh_.cell_halffaces(c)) {
                    std::vector<Eigen::Vector3d> triangle;
                    for (auto chv : mesh_.halfface_vertices(ch)) {
                        triangle.push_back(position(chv));
                    }
                    for (int i = 0; i < 3; i++) {
                        triangle[i] = scale_ * triangle[i] + (1 - scale_) * centroid;
                    }
                    Eigen::Vector3d n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]).normalized();
                    double nTd = n.dot(d);
                    if (nTd == 0) {
                        continue;
                    }
                    double lambda = (n.dot(triangle[0]) - n.dot(o)) / nTd;
                    if (lambda < 0) {
                        continue;
                    }
                    if (lambda < min) {
                        Eigen::Vector3d s = o + lambda * d;
                        for (int j = 0; j < 3; j++) {
                            triangle[j] -= s;
                        }
                        std::vector<Eigen::Vector3d> normals(3);
                        for (int j = 0; j < 3; j++) {
                            normals[j] = triangle[j].cross(triangle[(j + 1) % 3]);
                        }
                        if (normals[0].dot(normals[1]) >= 0 && normals[0].dot(normals[2]) >= 0) {
                            min = lambda;
                            type = 3;
                            idx = c.idx();
                        }
                    }
                }
            }
        }
#endif
        if (show_faces_ || twod_) {
            for (auto f : mesh_.faces()) {
                std::vector<Eigen::Vector3d> triangle;
                Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
                for (auto fv : mesh_.face_vertices(f)) {
                    triangle.push_back(position(fv));
                    centroid += position(fv);
                }
                centroid /= 3;
                for (int i = 0; i < 3; i++) {
                    triangle[i] = scale_ * triangle[i] + (1 - scale_) * centroid;
                }
                Eigen::Vector3d n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]).normalized();
                double nTd = n.dot(d);
                if (nTd == 0) {
                    continue;
                }
                double lambda = (n.dot(triangle[0]) - n.dot(o)) / nTd;
                if (lambda < 0) {
                    continue;
                }
                if (lambda < min) {
                    Eigen::Vector3d s = o + lambda * d;
                    for (int j = 0; j < 3; j++) {
                        triangle[j] -= s;
                    }
                    std::vector<Eigen::Vector3d> normals(3);
                    for (int j = 0; j < 3; j++) {
                        normals[j] = triangle[j].cross(triangle[(j + 1) % 3]);
                    }
                    if (normals[0].dot(normals[1]) >= 0 && normals[0].dot(normals[2]) >= 0) {
                        min = lambda;
                        type = 2;
                        idx = f.idx();
                    }
                }
            }
        }
        if (show_edges_ || twod_) {
            for (auto e : mesh_.edges()) {
                Halfedge h = mesh_.halfedge_handle(e, 0);
                Eigen::Vector3d a = position(mesh_.from_vertex_handle(h));
                Eigen::Vector3d b = position(mesh_.to_vertex_handle(h));
                double r = size_ / 1024;
                Eigen::Vector3d z = b - a;
                double length = z.norm();
                z /= length;
                Eigen::Vector3d y = z.cross(std::abs(z(0)) > 0.999 ? Eigen::Vector3d(0, 1, 0) : Eigen::Vector3d(1, 0, 0)).normalized();
                Eigen::Vector3d x = y.cross(z);
                Eigen::Matrix4d cw_matrix = Eigen::Matrix4d::Identity();
                cw_matrix.block<3, 1>(0, 0) = x;
                cw_matrix.block<3, 1>(0, 1) = y;
                cw_matrix.block<3, 1>(0, 2) = z;
                cw_matrix.block<3, 1>(0, 3) = a;
                Eigen::Matrix4d wc_matrix = cw_matrix.inverse();
                Eigen::Vector3d ca = (wc_matrix * o.homogeneous()).hnormalized();
                Eigen::Vector4d d_homogeneous;
                d_homogeneous << d, 0;
                Eigen::Vector3d cd = (wc_matrix * d_homogeneous).head(3);
                Eigen::Vector2d ca2 = ca.head(2);
                Eigen::Vector2d cd2 = cd.head(2);
                double cd2n = cd2.norm();
                if (cd2n == 0) {
                    continue;
                }
                cd2 /= cd2n;
                double aTd = ca2.dot(cd2);
                double radicant = aTd * aTd - ca2.dot(ca2) + r * r;
                if (radicant >= 0) {
                    double root = std::sqrt(radicant);
                    double lambda = (-aTd - root) / cd2n;
                    Eigen::Vector3d cs = ca + lambda * cd;
                    if (lambda < 0 || cs(2) < 0 || cs(2) > length) {
                        lambda = (-aTd + root) / cd2n;
                        cs = ca + lambda * cd;
                    }
                    if (lambda < 0 || cs(2) < 0 || cs(2) > length) {
                        continue;
                    }
                    if (lambda < min) {
                        min = lambda;
                        type = 1;
                        idx = e.idx();
                    }
                }
            }
        }
        if (show_vertices_ || twod_) {
            for (auto v : mesh_.vertices()) {
                Eigen::Vector3d c = position(v);
                double r = size_ / 512;
                Eigen::Vector3d a = o - c;
                double aTd = a.dot(d);
                double radicant = aTd * aTd - a.dot(a) + r * r;
                if (radicant >= 0) {
                    double root = std::sqrt(radicant);
                    double lambda = -aTd - root;
                    if (lambda < 0) {
                        lambda = -aTd + root;
                    }
                    if (lambda < 0) {
                        continue;
                    }
                    if (lambda < min) {
                        min = lambda;
                        type = 0;
                        idx = v.idx();
                    }
                }
            }
        }
        // Select
        if (idx != -1) {
            Eigen::Vector3d s = o + min * d;
            if (ImGui::IsKeyDown(ImGuiKey_LeftAlt)) {
#ifdef TET
                if (type == 3) {
                    auto sel = mesh_.property<Cell, bool>("sel", false);
                    Cell c(idx);
                    sel[c] = !sel[c];
                    if (sel[c]) {
                        std::cout << "Cell " << idx << std::endl;
                    }
                }
#endif
                if (type == 2) {
                    auto sel = mesh_.property<Face, bool>("sel", false);
                    Face f(idx);
                    sel[f] = !sel[f];
                    if (sel[f]) {
                        std::cout << "Face " << idx << std::endl;
                    }
                }
                if (type == 1) {
                    auto sel = mesh_.property<Edge, bool>("sel", false);
                    Edge e(idx);
                    sel[e] = !sel[e];
                    if (sel[e]) {
                        std::cout << "Edge " << idx << std::endl;
                    }
                }
                if (type == 0) {
                    auto sel = mesh_.property<Vertex, bool>("sel", false);
                    Vertex v(idx);
                    sel[v] = !sel[v];
                    if (sel[v]) {
                        std::cout << "Vertex " << idx << std::endl;
                    }
                }
                update();
            } else {
                s(2) = 0;
            }
            pick_callback_(type, idx, s);
        }
    }
    if (!twod_ && ImGui::IsWindowFocused() && ImGui::IsMouseDown(ImGuiMouseButton_Left) && !trans_ && !rot_cam_) {
        if (rot_) {
            Eigen::Vector2d m1 = to_screen(ImGui::GetMousePos(), size);
            if (m1 != m0_) {
                double r = (size.x > size.y ? size.y : size.x) / 2;
                Eigen::Vector3d v0;
                v0 << m0_, r;
                v0.normalize();
                Eigen::Vector3d v1;
                v1 << m1, r;
                v1.normalize();
                Eigen::Vector3d x = v0.cross(v1);
                double sin = x.norm();
                double cos = v0.dot(v1);
                Eigen::Vector3d n = x / sin;
                Eigen::Matrix3d xn;
                xn << 0, -n(2), n(1),
                    n(2), 0, -n(0),
                    -n(1), n(0), 0;
                Eigen::Matrix4d rotation_matrix = Eigen::Matrix4d::Identity();
                rotation_matrix.block<3, 3>(0, 0) = cos * Eigen::Matrix3d::Identity() + (1 - cos) * n * n.transpose() + sin * xn;
                transformation_matrix_ = rotation_matrix * transformation_matrix_;
                m0_ = m1;
            }
        } else {
            m0_ = to_screen(ImGui::GetMousePos(), size);
            rot_ = true;
        }
    } else {
        rot_ = false;
    }
    if (ImGui::IsMouseDown(ImGuiMouseButton_Middle) && !rot_ && !rot_cam_) {
        if (trans_) {
            Eigen::Vector2d m1 = to_screen(ImGui::GetMousePos(), size);
            if (m1 != m0_) {
                Eigen::Vector3d t;
                t << (m1 - m0_) * 2 * size_ / (size.x > size.y ? size.y : size.x) / fov_factor, 0;
                Eigen::Matrix4d translation_matrix = Eigen::Matrix4d::Identity();
                translation_matrix.block<3, 1>(0, 3) = t;
                transformation_matrix_ = translation_matrix * transformation_matrix_;
                m0_ = m1;
            }
        } else {
            m0_ = to_screen(ImGui::GetMousePos(), size);
            trans_ = true;
        }
    } else {
        trans_ = false;
    }
    if (ImGui::IsWindowHovered() && io.MouseWheel != 0) {
        if (ImGui::IsKeyDown(
#ifdef __APPLE__
                ImGuiKey_LeftSuper
#else
                ImGuiKey_LeftCtrl
#endif
                )) {
            if (io.MouseWheel > 0) {
                scale_ *= 1.148698355; // = 2^(1/5) -> 5 scrolls from 1 to 0.5
                if (scale_ > 1) {
                    scale_ = 1;
                }
            } else {
                scale_ /= 1.148698355;
            }
        } else {
            Eigen::Matrix4d scaling_matrix = Eigen::Matrix4d::Identity();
            if (io.MouseWheel > 0) {
                scaling_matrix.block<3, 3>(0, 0) *= 1.1;
            } else {
                scaling_matrix.block<3, 3>(0, 0) /= 1.1;
            }
            transformation_matrix_ = scaling_matrix * transformation_matrix_;
        }
    }
    if (!twod_ && ImGui::IsMouseDown(ImGuiMouseButton_Right) && !rot_ && !trans_) {
        if (rot_cam_) {
            Eigen::Vector2d m1 = to_screen(ImGui::GetMousePos(), size);
            if (m1 != m0_) {
                Eigen::Vector3d v0 = {0, 0, 1};
                Eigen::Vector2d delta = m1 - m0_;
                Eigen::Vector3d v1;
                v1 << delta, (size.x > size.y ? size.y : size.x) / 2;
                v1.normalize();
                Eigen::Vector3d x = v0.cross(v1);
                double sin = x.norm();
                double cos = v0.dot(v1);
                Eigen::Vector3d n = x / sin;
                Eigen::Matrix3d xn;
                xn << 0, -n(2), n(1),
                    n(2), 0, -n(0),
                    -n(1), n(0), 0;
                Eigen::Matrix4d rotation_matrix = Eigen::Matrix4d::Identity();
                rotation_matrix.block<3, 3>(0, 0) = cos * Eigen::Matrix3d::Identity() + (1 - cos) * n * n.transpose() + sin * xn;
                transformation_matrix_ = look_at_matrix.inverse() * rotation_matrix * look_at_matrix * transformation_matrix_;
                m0_ = m1;
            }
        } else {
            m0_ = to_screen(ImGui::GetMousePos(), size);
            rot_cam_ = true;
            glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        }
    } else {
        rot_cam_ = false;
        glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
    Eigen::Matrix4d translation_matrix = Eigen::Matrix4d::Identity();
    double distance = size_ * time;
    if (!twod_ && ImGui::IsKeyDown(ImGuiKey_W)) {
        translation_matrix(2, 3) += distance;
    }
    if (ImGui::IsKeyDown(ImGuiKey_A)) {
        translation_matrix(0, 3) += distance;
    }
    if (!twod_ && ImGui::IsKeyDown(ImGuiKey_S)) {
        translation_matrix(2, 3) -= distance;
    }
    if (ImGui::IsKeyDown(ImGuiKey_D)) {
        translation_matrix(0, 3) -= distance;
    }
    if (ImGui::IsKeyDown(ImGuiKey_Space) || twod_ && ImGui::IsKeyDown(ImGuiKey_W)) {
        translation_matrix(1, 3) -= distance;
    }
    if ((paper_ && ImGui::IsKeyDown(ImGuiKey_GraveAccent)) || (!paper_ && ImGui::IsKeyDown(ImGuiKey_LeftShift)) || twod_ && ImGui::IsKeyDown(ImGuiKey_S)) {
        translation_matrix(1, 3) += distance;
    }
    transformation_matrix_ = translation_matrix * transformation_matrix_;
    if (rotate_) {
        Eigen::Matrix4d rotation_matrix = Eigen::Matrix4d::Identity();
        double angle = 2 * M_PI * time / 4;
        double sin = std::sin(angle);
        double cos = std::cos(angle);
        rotation_matrix << cos, 0, sin, 0,
            0, 1, 0, 0,
            -sin, 0, cos, 0,
            0, 0, 0, 1;
        transformation_matrix_ = rotation_matrix * transformation_matrix_;
    }

    // Draw
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffer_);
    glBindTexture(GL_TEXTURE_2D, colorbuffer_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, pixel_ratio_ * size.x, pixel_ratio_ * size.y, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, colorbuffer_, 0);
    glBindRenderbuffer(GL_RENDERBUFFER, depthbuffer_);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, pixel_ratio_ * size.x, pixel_ratio_ * size.y);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, depthbuffer_);

    Eigen::Matrix4d model_view_matrix = look_at_matrix * transformation_matrix_;
    Eigen::Matrix4d normal_matrix = model_view_matrix.inverse().transpose();
    double near = size_ / 100;
    double far = size_ * 100;
    Eigen::Matrix4d projection_matrix = Eigen::Matrix4d::Identity() * fov_factor;
    if (size.x > size.y) {
        projection_matrix(0, 0) *= size.y / size.x;
    } else {
        projection_matrix(1, 1) *= size.x / size.y;
    }
    if (twod_) {
        projection_matrix /= size_;
        projection_matrix.block<2, 2>(2, 2) << 2 / (near - far), (near + far) / (near - far), 0, 1;
    } else {
        projection_matrix.block<2, 2>(2, 2) << (-far - near) / (far - near), -2 * far * near / (far - near), -1, 0;
    }
    float mvm[16];
    to_float_mat(model_view_matrix, mvm);
    glUniformMatrix4fv(glGetUniformLocation(shader_program_, "uModelViewMatrix"), 1, GL_FALSE, mvm);
    float nm[16];
    to_float_mat(normal_matrix, nm);
    glUniformMatrix4fv(glGetUniformLocation(shader_program_, "uNormalMatrix"), 1, GL_FALSE, nm);
    float pm[16];
    to_float_mat(projection_matrix, pm);
    glUniformMatrix4fv(glGetUniformLocation(shader_program_, "uProjectionMatrix"), 1, GL_FALSE, pm);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glUniform1i(glGetUniformLocation(shader_program_, "uTexture"), 0);
    glUniform1i(glGetUniformLocation(shader_program_, "uLight"), !twod_);
    glUniform1f(glGetUniformLocation(shader_program_, "uScale"), scale_);

    glViewport(0, 0, pixel_ratio_ * size.x, pixel_ratio_ * size.y);
    if (paper_) {
        glClearColor(1, 1, 1, 1);
    } else {
        glClearColor(0.1, 0.1, 0.1, 1);
    }
    glClearDepth(1);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

#ifdef TET
    if (show_cells_) {
        glBindVertexArray(cell_vao_);
        glDepthRange(0.001, 1);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell_index_buffer_);
        glDrawElements(GL_TRIANGLES, n_cell_indices_, GL_UNSIGNED_INT, 0);
        glDepthRange(0, 1);
    }
#endif
    if (show_faces_) {
        glBindVertexArray(face_vao_);
        if (use_texture_) {
            glUniform1i(glGetUniformLocation(shader_program_, "uUseTexture"), true);
        }
        glDepthRange(0.001, 1);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, face_index_buffer_);
        glDrawElements(GL_TRIANGLES, n_face_indices_, GL_UNSIGNED_INT, 0);
        glDepthRange(0, 1);
        glUniform1i(glGetUniformLocation(shader_program_, "uUseTexture"), false);
    }
    if (show_edges_) {
        glBindVertexArray(edge_vao_);
        glUniform1i(glGetUniformLocation(shader_program_, "uLight"), false);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_index_buffer_);
        glDrawElements(GL_LINES, n_edge_indices_, GL_UNSIGNED_INT, 0);
        glUniform1i(glGetUniformLocation(shader_program_, "uLight"), !twod_);
    }
    if (show_vertices_) {
        glBindVertexArray(vertex_instance_vao_);
        glDrawArraysInstanced(GL_TRIANGLES, 0, n_vertex_instance_indices_, n_vertex_instances_);
    }
    glBindVertexArray(extra_vao_);
    glDepthRange(0.001, 1);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, extra_index_buffer_);
    glDrawElements(GL_TRIANGLES, n_extra_indices_, GL_UNSIGNED_INT, 0);
    glDepthRange(0, 1);
    glUniform1i(glGetUniformLocation(shader_program_, "uLight"), false);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, extra_line_index_buffer_);
    glDrawElements(GL_LINES, n_extra_line_indices_, GL_UNSIGNED_INT, 0);
    glUniform1i(glGetUniformLocation(shader_program_, "uLight"), !twod_);

    ImGui::GetWindowDrawList()->AddImage(reinterpret_cast<ImTextureID>(colorbuffer_), ImVec2(0, 0), ImVec2(size.x, size.y), ImVec2(0, 1), ImVec2(1, 0));
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    ImGui::End();
    ImGui::PopStyleVar();

    prev_frame_ = curr_frame;
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window_);
}

Eigen::Vector3d Viewer::position(Vertex v) {
    if (parameter_space_) {
#ifndef TET
        auto tex = mesh_.property<Vertex, Eigen::Vector2d>("tex", {0.5, 0.5});
        Eigen::Vector3d p;
        p << tex[v], 0;
        return p;
#else
        auto parameter = mesh_.property<Vertex, Eigen::Vector3d>("parameter");
        return parameter[v];
#endif
    }
    return mesh_.position(v);
}

void Viewer::to_float_mat(const Eigen::Matrix4d& m, float array[16]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            array[4 * j + i] = m(i, j);
        }
    }
}

Eigen::Vector2d Viewer::to_screen(ImVec2 mouse_pos, ImVec2 size) {
    return Eigen::Vector2d(mouse_pos.x - size.x / 2, size.y / 2 - mouse_pos.y);
}
