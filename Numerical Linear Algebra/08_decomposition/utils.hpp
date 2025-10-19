#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <iostream>
#include <algorithm>

// from https://github.com/nothings/stb/tree/master
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace Eigen;
using DoubleMatrix = Matrix<double, Dynamic, Dynamic, RowMajor>;
using ImageMatrix = Matrix<unsigned char, Dynamic, Dynamic, RowMajor>;

DoubleMatrix clamp_matrix(const DoubleMatrix &m) {
    return m.cwiseMax(0.0).cwiseMin(255.0);
}

DoubleMatrix load_image_to_matrix(const char* input_image_path) {
    // Load the image using stb_image
    int width, height, channels;
    // for greyscale images force to load only one channel
    unsigned char* image_data = stbi_load(input_image_path, &width, &height, &channels, 1);
    if (!image_data) {
        std::cerr << "Error: Could not load image " << input_image_path << std::endl;
        exit(1);
    }

    std::cout << "Image " << input_image_path << " loaded." << std::endl;
    std::cout << "    Height:   " << height << std::endl;
    std::cout << "    Width:    " << width << std::endl;
    std::cout << "    Channels: " << channels << std::endl;
    
    // Prepare Eigen matrix
    DoubleMatrix image_matrix_xd(height, width);

    // Fill the matrix with image data
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int index = (i * width + j) * channels;  // 1 channel (Greyscale), 3 channels (RGB)
            image_matrix_xd(i, j) = static_cast<double>(image_data[index]);
        }
    }
    // Free memory!!!
    stbi_image_free(image_data);

    return image_matrix_xd;
}

void save_image(const DoubleMatrix &double_matrix, const std::string &output_image_path) {
    ImageMatrix image_matrix = clamp_matrix(double_matrix).cast<unsigned char>();
    
    int width = image_matrix.cols();
    int height = image_matrix.rows();

    // Save the image using stb_image_write
    if (stbi_write_png(output_image_path.c_str(), width, height, 1,
                       image_matrix.data(), width) == 0) {
        std::cerr << "Error: Could not save output image" << std::endl;
        exit(1);
    }

    std::cout << "Image saved to " << output_image_path << std::endl;
}

// This function applies filter by simple for loops (that is a bad way!).
// Only used to check the result.
DoubleMatrix apply_filter(const DoubleMatrix &image, const Matrix3d &filter) {
    // Applying filter in a non-effective way (using the for loops and formula at the top of challenge description)
    int height = image.rows();
    int width = image.cols();

    DoubleMatrix result = DoubleMatrix::Zero(height, width);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            double sum = 0.;

            // Apply filter to cell (i, j)
            for (int ki = 0; ki < 3; ki++) {
                for (int kj = 0; kj < 3; kj++) {
                    int ni = i + ki - 1;
                    int nj = j + kj - 1;
                    if (ni < 0 || ni >= height || nj < 0 || nj >= width) 
                        continue;
                    sum += filter(ki, kj) * image(ni, nj);
                }
            }
            
            result(i, j) = sum;
        }
    }

    // return clamp_matrix(result);
    return result;
}

SparseMatrix<double> construct_filter_matrix(int height, int width, const Matrix3d &filter) {
    int total_elements = height * width;
    SparseMatrix<double, RowMajor> filter_matrix(total_elements, total_elements);

    std::vector<Triplet<double>> triplet_list;

    for (int i = 0; i < total_elements; i++) {
        // the resulting image coordinates (row and col) of this index in filter matrix
        int result_row_id = i / width;
        int result_col_id = i % width;

        // set corresponding filter values
        for (int ki = 0; ki < 3; ki++) {
            for (int kj = 0; kj < 3; kj++) {
                int ni = result_row_id + ki - 1;
                int nj = result_col_id + kj - 1;
                if (ni < 0 || ni >= height || nj < 0 || nj >= width) 
                        continue;

                int cell_index = ni * width + nj;
            
                triplet_list.push_back(Triplet<double>(i, cell_index, filter(ki, kj)));
            }
        }
    }

    filter_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
    filter_matrix.makeCompressed();

    return filter_matrix;
}

VectorXd vector_from_matrix(const DoubleMatrix &matrix) {
    int height = matrix.rows();
    int width = matrix.cols();
    VectorXd v(height * width);

    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++){
            v(i * width + j) = matrix(i, j);
        }
    }

    return v;
}

bool is_matrix_symmetric(const SparseMatrix<double> &mat, double tolerance = 1e-10) {
    if (mat.rows() != mat.cols()) {
        return false;  // Must be square
    }
    
    return mat.isApprox(mat.transpose(), tolerance);
}

template<typename T>
void export_matrix(const std::string &filename, const SparseMatrix<T> &matrix) {
    saveMarket(matrix, filename);
    std::cout << "Matrix saved to " << filename << std::endl;
}

void export_vectorXd(const std::string &filename, const VectorXd &vec) {
    // Export vector in .mtx format
    int n = vec.size();
    // Eigen::saveMarketVector(b, "./rhs.mtx");
    FILE* out = fopen(filename.c_str(), "w");
    fprintf(out,"%%%%MatrixMarket vector coordinate real general\n");
    fprintf(out,"%d\n", n);
    for (int i=0; i<n; i++) {
        fprintf(out,"%d %f\n", i, vec(i));
    }
    fclose(out);
    
    std::cout << "Vector saved to " << filename << std::endl;
}

VectorXd load_vectorXd(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string line;
    
    // Skip comments and banner line
    while (std::getline(file, line)) {
        if (line[0] != '%') break;
    }
    
    // Read size (could be "n" or "n 1")
    int size, cols = 1;
    std::istringstream iss(line);
    iss >> size;
    if (iss >> cols) {
        // If there's a second number, it's a matrix format
        // Size is already correct
    }
    
    // Allocate vector
    Eigen::VectorXd vec(size);
    
    // Read values
    for (int i = 0; i < size; i++) {
        int dummyVal; file >> dummyVal;
        if (!(file >> vec(i))) {
            throw std::runtime_error("Error reading vector element " + std::to_string(i));
        }
    }
    
    file.close();

    return vec;
}

DoubleMatrix matrix_from_vector(const VectorXd &v, int height, int width) {
    DoubleMatrix matrix(height, width);

    // Fill the matrix with image data
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int index = (i * width + j);  // 1 channel (Greyscale)
            matrix(i, j) = v(index);
        }
    }

    return matrix;
}

VectorXd clamp_values(const VectorXd &v, double min_val, double max_val) {
    VectorXd result(v.rows());
    for (int i = 0; i < v.rows(); i++) {
        result(i) = std::clamp(v(i), min_val, max_val);
    }
    return result;
}

void print_task_start(int task_id) {
    std::cout << std::endl;
    std::cout << "******************" << std::endl;
    std::cout << "*     TASK " << task_id << (task_id < 10 ? " " : "") << "    *" << std::endl;
    std::cout << "******************" << std::endl;
}

std::string bool_str(bool b) {
    return (b ? "true" : "false");
}