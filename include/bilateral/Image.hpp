#ifndef BILATERALFILTER_IMAGE_HPP
#define BILATERALFILTER_IMAGE_HPP

#include <cmath>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <functional>

#include <glm/glm.hpp>

#include <bilateral/lodepng.hpp>


namespace bilateral {
    
    enum {
        NAIVE,
        SEPARABLE_KERNEL,
        BILATERAL_GRID
    };
    
    
    
    /**
     * Represents an Image with Pixel stored in a linear vector.
     *
     * @tparam N Number of channel per pixels.
     */
    template<int32_t N>
    class Image {
            static_assert((N == 1 || N == 3 || N == 4) && "Number of channel must be 1, 3 or 4");
        
        public:
            typedef glm::vec<N, uint8_t> Pixel;
        
        private:
            std::vector<Pixel> pixels; /**< Pixels array. */
            uint32_t height = 0; /**< Height of the Image. */
            uint32_t width = 0;  /**< Width of the Image. */
            uint32_t size = 0;   /**< Size of the Image. */
        
        public:
            
            ////////////////////////////////////////////////////////////////////
            ///////////////////////// CONSTRUCTORS /////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            /**
             * Initialize an empty and unusable Image.
             */
            Image() = default;
            
            
            /**
             * Load a PNG.
             *
             * @param filename Path to the PNG.
             */
            explicit Image(const char *filename);
            
            /**
             * Load a PNG.
             *
             * @param filename Path to the PNG.
             */
            explicit Image(const std::string &filename);
            
            
            /**
             * Initialize an empty image of size width * height.
             *
             * @param width Width of the image.
             * @param height Height of the image.
             */
            Image(uint32_t width, uint32_t height);
            
            /**
             * Initialize an image of size width * height with the given pixels.
             *
             * The array of pixel size must be equal to width * height.
             *
             * @param width Width if the image.g
             * @param height Height of the image.
             * @param pixels Pixels composing the image.
             */
            Image(uint32_t width, uint32_t height, const std::vector<uint8_t> &pixels);
            
            /**
             * Initialize an image of size width * height with the given pixels.
             *
             * The array of pixel size must be equal to width * height.
             *
             * @param width Width if the image.
             * @param height Height of the image.
             * @param pixels Pixels composing the image.
             */
            Image(uint32_t width, uint32_t height, const uint8_t *pixels);
            
            
            ////////////////////////////////////////////////////////////////////
            ///////////////////////// DATA ACCESS //////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            /**
             * Return a reference to a pixel.
             *
             * @tparam U Integral type use to as index of the pixel.
             *
             * @param index Index of the pixel.
             *
             * @return A reference to the pixel corresponding to index.
             */
            template<typename U>
            [[nodiscard]] Pixel &operator[](U index);
            
            /**
             * Return the value of a pixel.
             *
             * @tparam U Integral type use to as index of the pixel.
             *
             * @param index Index of the pixel.
             *
             * @return The value of the pixel corresponding to index.
             */
            template<typename U>
            [[nodiscard]] Pixel at(U index) const;
            
            /**
             * Return the vector of pixels composing the image.
             *
             * @return The vector of pixels composing the image.
             */
            [[nodiscard]] std::vector<uint8_t> getPixels() const;
            
            /**
             * Return the width of the image.
             *
             * @return The width of the image.
             */
            [[nodiscard]] uint32_t getWidth() const;
            
            /**
             * Return the height of the image.
             *
             * @return The height of the image
             */
            [[nodiscard]] uint32_t getHeight() const;
            
            /**
             * Return the number of pixel composing the image.
             *
             * @return The size of the image.
             */
            [[nodiscard]] uint32_t getSize() const;
            
            
            ////////////////////////////////////////////////////////////////////
            //////////////////////////// OTHERS ////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            /**
             * Create a copy of this image in dst, adding a border of size around the image.
             *
             * Border will be mirror reflection of the border elements (corresponding to BORDER_REFLECT on openCV).
             * Border's size cannot be greater the width or height.
             *
             * @param dst Destination of the bordered copy.
             * @param size Size of the border.
             */
            void borderedCopy(Image<N> &dst, uint32_t size) const;
            
            /**
             * Save the image to path as a PNG.
             *
             * @param path Path to the output image.
             */
            void save(const char *path) const;
            
            /**
             * Compute a 1D gaussian kernel of diameter radius * 2 + 1.
             *
             * @param kernel Destination of the computed kernel.
             * @param sigma Sigma used to compute the kernel.
             * @param radius Radius of the kernel.
             */
            static void gaussianKernel1D(std::vector<glm::vec<N, double>> &kernel, double sigma, int32_t radius);
            
            /**
             * Compute a 2D gaussian kernel of diameter radius * 2 + 1.
             *
             * @param kernel Destination of the computed kernel.
             * @param sigma Sigma used to compute the kernel.
             * @param radius Radius of the kernel.
             */
            static void gaussianKernel2D(std::vector<glm::vec<N, double>> &kernel, double sigma, int32_t radius);
            
            /**
             * Compute a 3D gaussian kernel of diameter radius * 2 + 1.
             *
             * @param kernel Destination of the computed kernel.
             * @param sSpace Sigma used for the space weight.
             * @param sRange Sigma used for the range weight.
             * @param radius Radius of the kernel..
             */
            static void gaussianKernel3D(std::vector<glm::vec<N, double>> &kernel, double sSpace, double sRange,
                                         int32_t radius);
            
            /**
             * Compute the gaussian(x) where x is the difference between p1 and p2.
             *
             * @param p1 First pixel.
             * @param p2 Second pixel.
             * @param gaussian Gaussian function.
             *
             * @return The computed gaussian value.
             */
            static glm::vec<N, double> computeRangeGaussian(const Pixel &p1, const Pixel &p2,
                                                            const std::function<double(double)> &gaussian);
            
            /**
             * Compute in dst the result of a bilateral filter on src using different sigma for the space and range factor.
             *
             * Brute force implementation of the bilateral filter.
             *
             * @param dst Destination of the result of the bilateral filter on src.
             * @param sSpace Sigma used for the space weight.
             * @param sRange Sigma used for the range weight.
             */
            void naive(Image<N> &dst, double sSpace, double sRange);
            
            
            /**
             * Compute in dst the result of a bilateral filter on src using different sigma for the space and range factor.
             *
             * Faster implementation but introduce some inaccuracy.
             *
             * @param dst Destination of the result of the bilateral filter on src.
             * @param sSpace Sigma used for the space weight.
             * @param sRange Sigma used for the range weight.
             */
            void separableKernel(Image<N> &dst, double sSpace, double sRange);
            
            /**
             * Compute in dst the result of a bilateral filter on src using different sigma for the space and range factor.
             *
             * Use a bilateral grid, as introduce by Jiawen Chen, Sylvain Paris, Frédo Durand in their paper,
             * "Real-time Edge-Aware Image Processing with the Bilateral Grid".
             *
             * @param dst Destination of the result of the bilateral filter on src.
             * @param sSpace Sigma used for the space weight.
             * @param sRange Sigma used for the range weight.
             * @param sampleSpace Sample rate used for the spatial axes.
             * @param sampleRange Sample rate used for the range axis.
             */
            void bilateralGrid(Image<N> &dst, double sSpace, double sRange, double sampleSpace, double sampleRange);
    };
    
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// TYPES ALIASES ///////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    using Image1u8 = Image<1>;
    
    using Image3u8 = Image<3>;
    
    using Image4u8 = Image<4>;
    
    
    ////////////////////////////////////////////////////////////////////////////
    /////////////////////////// TEMPLATE DEFINITION ////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    template<int32_t N>
    Image<N>::Image(const std::string &filename) :
            Image(filename.c_str()) {
    }
    
    
    template<int32_t N>
    Image<N>::Image(uint32_t width, uint32_t height) :
            height(height), width(width), size(height * width) {
        this->pixels.resize(this->size);
    }
    
    
    template<int32_t N>
    Image<N>::Image(uint32_t width, uint32_t height, const std::vector<uint8_t> &pixels) :
            pixels(pixels), height(height), width(width), size(height * width) {
        if (this->pixels.size() != this->size) {
            std::runtime_error("Error: The size of the given pixel array is not equal to width * height");
        }
    }
    
    
    template<int32_t N>
    Image<N>::Image(uint32_t width, uint32_t height, const uint8_t *pixels) :
            pixels(pixels, pixels + width * height), height(height), width(width), size(height * width) {
        if (this->pixels.size() != this->size) {
            std::runtime_error("Error: The size of the given pixel array is not equal to width * height");
        }
    }
    
    
    template<int32_t N>
    template<typename U>
    auto Image<N>::operator[](U index) -> Image<N>::Pixel & {
        static_assert(std::is_integral<U>::value, "Integral type required.");
        return this->pixels[static_cast<uint32_t >(index)];
    }
    
    
    template<int32_t N>
    template<typename U>
    auto Image<N>::at(U index) const -> Image<N>::Pixel {
        static_assert(std::is_integral<U>::value, "Integral type required.");
        return this->pixels.at(static_cast<uint32_t >(index));
    }
    
    
    template<int32_t N>
    std::vector<uint8_t> Image<N>::getPixels() const {
        return this->pixels;
    }
    
    
    template<int32_t N>
    uint32_t Image<N>::getWidth() const {
        return this->width;
    }
    
    
    template<int32_t N>
    uint32_t Image<N>::getHeight() const {
        return this->height;
    }
    
    
    template<int32_t N>
    uint32_t Image<N>::getSize() const {
        return this->size;
    }
    
    
    template<int32_t N>
    void Image<N>::borderedCopy(Image<N> &dst, uint32_t size) const {
        assert(size < width && size < height && "The size of the border cannot be greater than width or height");
        
        uint32_t width = this->width + size + size;
        uint32_t height = this->height + size + size;
        dst = Image<N>(width, height);
        uint32_t dstY, srcY;
        
        for (uint32_t y = 0; y < this->height; y++) {
            dstY = (y + size) * width;
            srcY = y * this->width;
            for (uint32_t x = 0; x < this->width; x++) {
                dst[dstY + x + size] = this->pixels[srcY + x];
            }
        }
    }
    
    //////////////////////////// GAUSSIAN KERNEL ///////////////////////////////
    
    template<int32_t N>
    void Image<N>::gaussianKernel1D(std::vector<glm::vec<N, double>> &kernel, double sigma, int32_t radius) {
        const double sigma2 = sigma * sigma;
        const double factor = 1. / (2. * M_PI * sigma2);
        const double divisor = 2. * sigma2;
        
        // Gaussian function with precomputed factor and divisor
        auto sGaussian = [factor, divisor](double value) {
            return factor * std::exp(-(value * value) / divisor);
        };
        
        const int32_t diameter = radius * 2 + 1;
        uint32_t index = 0;
        glm::vec<N, double> sum = glm::vec<N, double>(0);
        
        // Initialize the kernel size
        kernel.resize(static_cast<uint32_t >(diameter));
        
        
        // Fill the kernel
        for (int32_t y = -radius; y <= radius; y++, index++) {
            kernel[index] = glm::vec<N, double>(sGaussian(std::abs(y)));
            sum += glm::vec<N, double>(kernel[index]);
        }
        // Normalize the probability mass outside the kernel evenly to all pixels within the kernel
        for (glm::vec<N, double> &x : kernel) {
            x /= sum;
        }
    }
    
    
    template<int32_t N>
    void Image<N>::gaussianKernel2D(std::vector<glm::vec<N, double>> &kernel, double sigma, int32_t radius) {
        const double sigma2 = sigma * sigma;
        const double factor = 1. / (2. * M_PI * sigma2);
        const double divisor = 2. * sigma2;
        
        // Gaussian function with precomputed factor and divisor
        auto sGaussian = [factor, divisor](double value) {
            return factor * std::exp(-(value * value) / divisor);
        };
        
        const int32_t diameter = radius * 2 + 1;
        uint32_t index = 0;
        glm::vec<N, double> sum = glm::vec<N, double>(0);
        
        // Initialize the kernel size
        kernel.resize(static_cast<uint32_t >(diameter * diameter));
        
        // Fill the kernel
        for (int32_t y = -radius; y <= radius; y++) {
            for (int32_t x = -radius; x <= radius; x++, index++) {
                kernel[index] = glm::vec<N, double>(sGaussian(std::sqrt(y * y + x * x)));
                sum += glm::vec<N, double>(kernel[index]);
            }
        }
        // Normalize the probability mass outside the kernel evenly to all pixels within the kernel
        for (glm::vec<N, double> &x : kernel) {
            x /= sum;
        }
    }
    
    
    template<int32_t N>
    void Image<N>::gaussianKernel3D(std::vector<glm::vec<N, double>> &kernel, double sSpace, double sRange,
                                    int32_t radius) {
        const double sSigma2 = sSpace * sSpace;
        const double rDivisor = 2. * sRange * sRange;
        const double sDivisor = 2. * sSigma2;
        const double factor = 1. / (2. * M_PI * sSigma2 * sRange);
        
        const auto gaussian = [factor, sDivisor, rDivisor](double space, double range) {
            return factor * std::exp(-((space * space) / sDivisor + range / rDivisor));
        };
        
        const int32_t diameter = radius * 2 + 1;
        glm::vec<N, double> sum = glm::vec<N, double>(0);
        uint32_t index = 0;
        
        // Initialize the kernel size
        kernel.resize(static_cast<uint32_t >(diameter * diameter * diameter));
        
        // Fill the kernel
        for (int32_t z = -radius; z <= radius; z++) {
            for (int32_t y = -radius; y <= radius; y++) {
                for (int32_t x = -radius; x <= radius; x++, index++) {
                    kernel[index] = glm::vec<N, double>(gaussian(std::sqrt(y * y + x * x)), std::abs(z));
                    sum += glm::vec<N, double>(kernel[index]);
                }
            }
        }
        
        // Normalize the probability mass outside the kernel evenly to all pixels within the kernel
        for (glm::vec<N, double> &x : kernel) {
            x /= sum;
        }
    }
    
    
    template<int32_t N>
    glm::vec<N, double> Image<N>::computeRangeGaussian(const Pixel &p1, const Pixel &p2,
                                                       const std::function<double(double)> &gaussian) {
        glm::vec<N, double> g;
        
        for (int32_t i = 0; i < N; i++) {
            g[i] = gaussian(std::abs(p1[0] - p2[0]));
        }
        
        return g;
    }
    
    ///////////////////////////////// NAIVE ////////////////////////////////////
    
    template<int32_t N>
    void Image<N>::naive(Image<N> &dst, double sSpace, double sRange) {
        assert(sRange > 0 && "Color sigma must be greater than 0");
        assert(sSpace > 0 && "Space sigma must be greater than 0");
        
        const int32_t radius = static_cast<const int32_t>(std::round(sSpace * 1.5));
        const int32_t diameter = radius * 2 + 1;
        Image<N> bordered;
        
        this->borderedCopy(bordered, static_cast<uint32_t>(radius));
        dst = Image(this->getWidth(), this->getHeight());
        
        std::vector<glm::vec<N, double>> sKernel;
        gaussianKernel2D(sKernel, sSpace, radius);
        
        const double sigma2 = sRange * sRange;
        const double factor = 1. / (2. * M_PI * sigma2);
        const double divisor = 2. * sigma2;
        std::function<double(double)> rGaussian = std::bind(
                [](double factor, double divisor, double value) {
                    return factor * std::exp(-(value * value) / divisor);
                },
                factor, divisor, std::placeholders::_1
        );
        
        int32_t widthb = static_cast<int32_t>(bordered.getWidth());
        int32_t width = static_cast<int32_t>(dst.getWidth());
        int32_t height = static_cast<int32_t>(dst.getHeight());
        glm::vec<N, double> sum, wp, kmul;
        uint32_t indexk, center, index;
        
        for (int32_t y = 0; y < height; y++) {
            for (int32_t x = 0; x < width; x++) {
                sum = glm::vec<N, double>(0);
                wp = glm::vec<N, double>(0);
                
                center = static_cast<uint32_t>((y + radius) * widthb + x + radius);
                for (int32_t ky = -radius; ky <= radius; ky++) {
                    for (int32_t kx = -radius; kx <= radius; kx++) {
                        indexk = static_cast<uint32_t>((ky + radius) * diameter + kx + radius);
                        index = static_cast<uint32_t>((y + ky + radius) * widthb + x + kx + radius);
                        kmul = sKernel[indexk] * computeRangeGaussian(bordered[index], bordered[center], rGaussian);
                        wp += kmul;
                        sum += kmul * glm::vec<N, double>(bordered[index]);
                    }
                }
                
                dst[y * width + x] = static_cast<Pixel>(1. / wp * sum);
            }
        }
    }
    
    
    
    /////////////////////////// SEPARABLE KERNEL ///////////////////////////////
    
    
    template<int32_t N>
    void Image<N>::separableKernel(Image<N> &dst, double sSpace, double sRange) {
        assert(sRange > 0 && "Color sigma must be greater than 0");
        assert(sSpace > 0 && "Space sigma must be greater than 0");
        
        const int32_t radius = static_cast<const int32_t>(std::round(sSpace * 1.5));
        Image<N> bordered;
        
        this->borderedCopy(bordered, static_cast<uint32_t>(radius));
        dst = Image(this->getWidth(), this->getHeight());
        
        std::vector<glm::vec<N, double>> sKernel;
        gaussianKernel1D(sKernel, sSpace, radius);
        
        const double sigma2 = sRange * sRange;
        const double factor = 1. / (2. * M_PI * sigma2);
        const double divisor = 2. * sigma2;
        std::function<double(double)> rGaussian = std::bind(
                [](double factor, double divisor, double value) {
                    return factor * std::exp(-(value * value) / divisor);
                },
                factor, divisor, std::placeholders::_1
        );
        
        int32_t widthb = static_cast<int32_t>(bordered.getWidth());
        int32_t width = static_cast<int32_t>(dst.getWidth());
        int32_t height = static_cast<int32_t>(dst.getHeight());
        glm::vec<N, double> sum, wp, kmul;
        uint32_t indexk, center, index;
        
        // Apply on row
        for (int32_t y = 0; y < height; y++) {
            for (int32_t x = 0; x < width; x++) {
                sum = glm::vec<N, double>(0);
                wp = glm::vec<N, double>(0);
                
                center = static_cast<uint32_t>((y + radius) * widthb + x + radius);
                for (int32_t kx = -radius; kx <= radius; kx++) {
                    indexk = static_cast<uint32_t>(kx + radius);
                    index = static_cast<uint32_t>((y + radius) * widthb + x + kx + radius);
                    kmul = sKernel[indexk] * computeRangeGaussian(bordered[index], bordered[center], rGaussian);
                    wp += kmul;
                    sum += kmul * glm::vec<N, double>(bordered[index]);
                }
                
                dst[y * width + x] = static_cast<Pixel>(1. / wp * sum);
            }
        }
        
        dst.borderedCopy(bordered, static_cast<uint32_t>(radius));
        dst = Image(this->getWidth(), this->getHeight());
        
        // Apply on column
        for (int32_t y = 0; y < height; y++) {
            for (int32_t x = 0; x < width; x++) {
                sum = glm::vec<N, double>(0);
                wp = glm::vec<N, double>(0);
                
                center = static_cast<uint32_t>((y + radius) * widthb + x + radius);
                for (int32_t ky = -radius; ky <= radius; ky++) {
                    indexk = static_cast<uint32_t>(ky + radius);
                    index = static_cast<uint32_t>((y + ky + radius) * widthb + x + radius);
                    kmul = sKernel[indexk] * computeRangeGaussian(bordered[index], bordered[center], rGaussian);
                    wp += kmul;
                    sum += kmul * glm::vec<N, double>(bordered[index]);
                }
                
                dst[y * width + x] = static_cast<Pixel>(1. / wp * sum);
            }
        }
    }
    
    //////////////////////////// BILATERAL GRID ////////////////////////////////
    
    template<int32_t N>
    static void borderedGridCopy(const std::vector<glm::vec<N + 1, double>> &src,
                                 std::vector<glm::vec<N + 1, double>> &dst, uint32_t xSize, uint32_t ySize,
                                 uint32_t zSize, uint32_t borderSize) {
        uint32_t width = xSize + borderSize + borderSize;
        uint32_t height = ySize + borderSize + borderSize;
        uint32_t depth = zSize + borderSize + borderSize;
        uint32_t dstY, srcY, dstZ, srcZ;
        
        dst.resize(width * height * depth);
        
        for (uint32_t z = 0; z < zSize; z++) {
            dstZ = (z + borderSize) * zSize;
            srcZ = z * zSize;
            
            for (uint32_t y = 0; y < ySize; y++) {
                dstY = (y + borderSize + dstZ) * ySize;
                srcY = (y + srcZ) * ySize;
                
                for (uint32_t x = 0; x < xSize; x++) {
                    dst[dstY + x + borderSize] = src[srcY + x];
                }
            }
        }
    }
    
    
    template<int32_t N>
    void Image<N>::bilateralGrid(Image<N> &dst, double sSpace, double sRange, double sampleSpace,
                                 double sampleRange) {
        assert(sSpace > 0 && "Space sigma must be greater than 0");
        assert(sRange > 0 && "Color sigma must be greater than 0");
        assert(sampleSpace > 0 && "Sample space must be greater than 0");
        assert(sampleRange > 0 && "Sample range must be greater than 0");
        
        auto grid = std::vector<glm::vec<N + 1, double>>(this->width * this->height * 255);
        uint32_t i, iGrid;
        
        for (uint32_t y = 0; y < this->height; y++) {
            for (uint32_t x = 0; x < this->width; x++) {
                i = x + y * this->width;
                iGrid = x + this->width * (y + 255 * this->pixels[i].x);
                
                grid[iGrid] = std::vector<glm::vec<N + 1, uint8_t>>(this->pixels[i], 1);
            }
        }
        
        const int32_t radius = static_cast<const int32_t>(std::round(sSpace * 1.5));
        
        Image<N> bordered;
        this->borderedGridCopy(grid, bordered, this->width, this->height, 255, radius);
        
        std::vector<glm::vec<N + 1, double>> sKernel;
        gaussianKernel3D(sKernel, sSpace, sRange, radius);
    }
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////// TEMPLATE SPECIALIZATION //////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    template<>
    Image<1>::Image(const char *filename) {
        std::vector<uint8_t> raw;
        
        decodePNG(raw, this->width, this->height, filename);
        
        this->size = this->width * this->height;
        this->pixels.resize(this->size);
        // lodepng loads in RGBA, we need to convert it in gray scale.
        for (uint32_t i = 0; i < this->size; i++) {
            this->pixels[i].x = (raw[i * 4] + raw[i * 4 + 1] + raw[i * 4 + 2]) / 3;
        }
    }
    
    
    template<>
    void Image<1>::save(const char *path) const {
        std::vector<uint8_t> raw(this->size * 4);
        
        for (uint32_t i = 0; i < this->size; i++) {
            raw[i * 4] = this->pixels[i].x;
            raw[i * 4 + 1] = this->pixels[i].x;
            raw[i * 4 + 2] = this->pixels[i].x;
            raw[i * 4 + 3] = 255;
        }
        
        encodePNG(path, raw, this->width, this->height);
    }
    
    
    template<>
    Image<3>::Image(const char *filename) {
        std::vector<uint8_t> raw;
        
        decodePNG(raw, this->width, this->height, filename);
        
        this->size = this->width * this->height;
        this->pixels.resize(this->size);
        // lodepng loads in RGBA, we need to convert it in RGBA.
        for (uint32_t i = 0; i < this->size; i++) {
            this->pixels[i].r = raw[i * 4];
            this->pixels[i].g = raw[i * 4 + 1];
            this->pixels[i].b = raw[i * 4 + 2];
        }
    }
    
    
    template<>
    void Image<3>::save(const char *path) const {
        std::vector<uint8_t> raw(this->size * 4);
        
        for (uint32_t i = 0; i < this->size; i++) {
            raw[i * 4] = this->pixels[i].r;
            raw[i * 4 + 1] = this->pixels[i].g;
            raw[i * 4 + 2] = this->pixels[i].b;
            raw[i * 4 + 3] = 255;
        }
        
        encodePNG(path, raw, this->width, this->height);
    }
    
    
    template<>
    Image<4>::Image(const char *filename) {
        std::vector<uint8_t> raw;
        
        decodePNG(raw, this->width, this->height, filename);
        
        this->size = this->width * this->height;
        this->pixels.resize(this->size);
        
        for (uint32_t i = 0; i < this->size; i++) {
            this->pixels[i].r = raw[i * 4];
            this->pixels[i].g = raw[i * 4 + 1];
            this->pixels[i].b = raw[i * 4 + 2];
            this->pixels[i].a = raw[i * 4 + 3];
        }
    }
    
    
    template<>
    void Image<4>::save(const char *path) const {
        std::vector<uint8_t> raw(this->size * 4);
        
        for (uint32_t i = 0; i < this->size; i++) {
            raw[i * 4] = this->pixels[i].r;
            raw[i * 4 + 1] = this->pixels[i].g;
            raw[i * 4 + 2] = this->pixels[i].b;
            raw[i * 4 + 3] = this->pixels[i].a;
        }
        
        encodePNG(path, raw, this->width, this->height);
    }
}

#endif // BILATERALFILTER_IMAGE_HPP
