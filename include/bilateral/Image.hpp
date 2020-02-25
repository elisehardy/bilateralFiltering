#ifndef BILATERALFILTER_IMAGE_HPP
#define BILATERALFILTER_IMAGE_HPP

#include <cmath>
#include <string>
#include <utility>
#include <vector>
#include <iostream>

#include <glm/glm.hpp>

#include <bilateral/lodepng.hpp>
#include <functional>


namespace bilateral {
    
    /**
     * Represents an Image with Pixel stored in a linear vector.
     *
     * @tparam N Number of channel per pixels.
     * @tparam T Type used to represent a channel.
     */
    template<int32_t N, typename T>
    class Image {
            static_assert((N == 1 || N == 3 || N == 4) && "Number of channel must be 1, 3 or 4");
            static_assert(std::is_arithmetic<T>::value, "Arithmetic type required.");
        
        public:
            typedef glm::vec<N, T> Pixel;
        
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
            Image(uint32_t width, uint32_t height, std::vector<T> pixels);
            
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
            [[nodiscard]] std::vector<glm::vec<N, T>> getPixels() const;
            
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
             * Compute a 2D gaussian kernel of diameter radius * 2 + 1.
             *
             * @param kernel Destination of the computed kernel.
             * @param sigma Sigma used to compute the kernel.
             * @param radius Radius of the kernel.
             */
            static void gaussianKernel2D(std::vector<glm::vec<N, double>> &kernel, double sigma, int32_t radius);
            
            /**
             * Compute a 1D gaussian kernel of diameter radius * 2 + 1.
             *
             * @param kernel Destination of the computed kernel.
             * @param sigma Sigma used to compute the kernel.
             * @param radius Radius of the kernel.
             */
            static void gaussianKernel1D(std::vector<double> &kernel, double sigma, int32_t radius);
            
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
             * Create a copy of this image in dst, adding a border of size around the image.
             *
             * Border will be mirror reflection of the border elements (corresponding to BORDER_REFLECT on openCV).
             * Border's size cannot be greater the width or height.
             *
             * @param dst Destination of the bordered copy.
             * @param size Size of the border.
             */
            void borderedCopy(Image<N, T> &dst, uint32_t size) const;
            
            /**
            * Compute in dst the result of a bilateral filter on src using different sigma for the space and range factor.
            *
            * @param dst Destination of the result of the bilateral filter on src.
            * @param sSpace Sigma use for the space weight.
            * @param sRange Sigma use for the range weight.
            */
            void naive(Image<N, T> &dst, double sSpace, double sRange);
            
            /**
             * Save the image to path as a PNG.
             *
             * @param path Path to the output image.
             */
            void save(const char *path) const;
    };
    
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// TYPES ALIASES ///////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    using Image1u8 = Image<1, uint8_t>;
    
    using Image3u8 = Image<3, uint8_t>;
    
    using Image1f32 = Image<1, float>;
    
    using Image3f32 = Image<3, float>;
    
    using Image1f64 = Image<1, double>;
    
    using Image3f64 = Image<1, double>;
    
    
    ////////////////////////////////////////////////////////////////////////////
    /////////////////////////// TEMPLATE DEFINITION ////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    template<int32_t N, typename T>
    Image<N, T>::Image(const std::string &filename) :
            Image(filename.c_str()) {
    }
    
    
    template<int32_t N, typename T>
    Image<N, T>::Image(uint32_t width, uint32_t height) :
            height(height), width(width), size(height * width) {
        this->pixels.resize(this->size);
    }
    
    
    template<int32_t N, typename T>
    Image<N, T>::Image(uint32_t width, uint32_t height, std::vector<T> pixels) :
            pixels(std::move(pixels)), height(height), width(width), size(height * width) {
        if (this->pixels.size() != this->size) {
            std::runtime_error("Error: The size of the given pixel array is not equal to width * height");
        }
    }
    
    
    template<int32_t N, typename T>
    Image<N, T>::Image(uint32_t width, uint32_t height, const uint8_t *pixels) :
            pixels(pixels, pixels + width * height), height(height), width(width), size(height * width) {
        if (this->pixels.size() != this->size) {
            std::runtime_error("Error: The size of the given pixel array is not equal to width * height");
        }
    }
    
    
    template<int32_t N, typename T>
    template<typename U>
    auto Image<N, T>::operator[](U index) -> Image<N, T>::Pixel & {
        static_assert(std::is_integral<U>::value, "Integral type required.");
        return this->pixels[static_cast<uint32_t >(index)];
    }
    
    
    template<int32_t N, typename T>
    template<typename U>
    auto Image<N, T>::at(U index) const -> Image<N, T>::Pixel {
        static_assert(std::is_integral<U>::value, "Integral type required.");
        return this->pixels.at(static_cast<uint32_t >(index));
    }
    
    
    template<int32_t N, typename T>
    std::vector<glm::vec<N, T>> Image<N, T>::getPixels() const {
        return this->pixels;
    }
    
    
    template<int32_t N, typename T>
    uint32_t Image<N, T>::getWidth() const {
        return this->width;
    }
    
    
    template<int32_t N, typename T>
    uint32_t Image<N, T>::getHeight() const {
        return this->height;
    }
    
    
    template<int32_t N, typename T>
    uint32_t Image<N, T>::getSize() const {
        return this->size;
    }
    
    
    template<int32_t N, typename T>
    void Image<N, T>::gaussianKernel1D(std::vector<double> &kernel, double sigma, int32_t radius) {
        const double sigma2 = sigma * sigma;
        const double factor = 1. / (2. * M_PI * sigma2);
        const double divisor = 2. * sigma2;
        
        // Gaussian function with precomputed factor and divisor
        auto sGaussian = [factor, divisor](double value) {
            return factor * std::exp(-(value * value) / divisor);
        };
        
        const int32_t diameter = radius * 2 + 1;
        uint32_t index = 0;
        double sum = 0;
        
        // Initialize the kernel size
        kernel.resize(static_cast<uint32_t >(diameter));
        
        // Fill the kernel with the gaussian value corresponding to the center of the kernel
        for (int32_t y = -radius; y <= radius; y++) {
            kernel[index] = sGaussian(std::abs(y));
            sum += kernel[index];
        }
        // Normalize the probability mass outside the kernel evenly to all pixels within the kernel
        for (double &x : kernel) {
            x /= sum;
        }
    }
    
    
    template<int32_t N, typename T>
    void Image<N, T>::gaussianKernel2D(std::vector<glm::vec<N, double>> &kernel, double sigma, int32_t radius) {
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
        
        // Fill the kernel with the gaussian value corresponding to the center of the kernel
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
    
    
    template<int32_t N, typename T>
    glm::vec<N, double> Image<N, T>::computeRangeGaussian(const Pixel &p1, const Pixel &p2,
                                                          const std::function<double(double)> &gaussian) {
        glm::vec<N, double> g;
        
        for (int32_t i = 0; i < N; i++) {
            g[i] = gaussian(std::abs(p1[0] - p2[0]));
        }
        
        return g;
    }
    
    
    template<int32_t N, typename T>
    void Image<N, T>::borderedCopy(Image<N, T> &dst, uint32_t size) const {
        assert(size < width && size < height && "The size of the border cannot be greater than width or height");
        
        uint32_t width = this->width + size + size;
        uint32_t height = this->height + size + size;
        uint32_t dstY, srcY;
        
        dst = Image<N, T>(width, height);
        
        // Center
        for (uint32_t y = 0; y < this->height; y++) {
            dstY = (y + size) * width;
            srcY = y * this->width;
            for (uint32_t x = 0; x < this->width; x++) {
                dst[dstY + x + size] = this->pixels[srcY + x];
            }
        }
        
        // Top border
        for (uint32_t y = 0; y < size; y++) {
            dstY = y * width;
            srcY = (size - y) * this->width;
            for (uint32_t x = size; x < this->width + size; x++) {
                dst[dstY + x] = this->pixels[srcY + x - size];
            }
        }
        
        // Bottom border
        for (uint32_t y = this->height + size, tmp = 0; y < height; y++, tmp++) {
            dstY = y * width;
            srcY = (this->height - tmp - 1) * this->width;
            for (uint32_t x = size; x < this->width + size; x++) {
                dst[dstY + x] = this->pixels[srcY + x - size];
            }
        }
        
        // Left border
        for (uint32_t y = 0; y < height; y++) {
            dstY = y * width;
            for (uint32_t x = 0, tmp = 0; x < size; x++, tmp++) {
                dst[dstY + x] = dst[dstY + size + size - tmp];
            }
        }
        
        // Right border
        for (uint32_t y = 0; y < height; y++) {
            dstY = y * width;
            for (uint32_t x = this->width + size, tmp = 0; x < width; x++, tmp++) {
                dst[dstY + x] = dst[dstY + this->width + size - tmp - 1];
            }
        }
    }
    
    
    template<int32_t N, typename T>
    void Image<N, T>::naive(Image<N, T> &dst, double sSpace, double sRange) {
        assert(sRange > 0 && "Color sigma must be greater than 0");
        assert(sSpace > 0 && "Space sigma must be greater than 0");
        
        const int32_t radius = static_cast<const int32_t>(std::round(sSpace * 1.5));
        const int32_t diameter = radius * 2 + 1;
        Image<N, T> bordered;
        
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
        uint32_t indexk, center, index;
        glm::vec<N, double> sum, wp, kmul;
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
    
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////// TEMPLATE SPECIALIZATION //////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    template<>
    Image<1, uint8_t>::Image(const char *filename) {
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
    void Image<1, uint8_t>::save(const char *path) const {
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
    Image<3, uint8_t>::Image(const char *filename) {
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
    void Image<3, uint8_t>::save(const char *path) const {
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
    Image<4, uint8_t>::Image(const char *filename) {
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
    void Image<4, uint8_t>::save(const char *path) const {
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
