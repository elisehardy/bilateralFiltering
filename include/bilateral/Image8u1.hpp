#ifndef BILATERALFILTER_IMAGE8U1_HPP
#define BILATERALFILTER_IMAGE8U1_HPP

#include <string>
#include <utility>
#include <vector>
#include <iostream>


namespace bilateral {
    
    /**
     * Represents a gray-scale image.
     *
     * Each pixel is linearly stored (row-major order) as unsigned char.
     */
    class Image8u1 {
        
        private:
            std::vector<uint8_t> pixels;
            uint32_t height = 0;
            uint32_t width = 0;
            uint32_t size = 0;
        
        public:
            
            ////////////////////////////////////////////////////////////////////
            ///////////////////////// CONSTRUCTORS /////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            /**
             * Initialize an empty and unusable Image.
             */
            Image8u1() = default;
            
            /**
             * Load a PNG as a gray-scale image.
             *
             * @param filename Path to the PNG.
             */
            explicit Image8u1(const char *filename);
            
            /**
             * Load a PNG as a gray-scale image.
             *
             * @param filename Path to the PNG.
             */
            explicit Image8u1(const std::string &filename);
            
            /**
             * Initialize an empty image of size width * height.
             *
             * @param width Width if the image.
             * @param height Height of the image.
             */
            Image8u1(uint32_t width, uint32_t height);
            
            /**
             * Initialize an image of size width * height with the given pixels.
             *
             * The array of pixel size must be equal to width * height.
             *
             * @param width Width if the image.g
             * @param height Height of the image.
             * @param pixels Pixels composing the image.
             */
            Image8u1(uint32_t width, uint32_t height, std::vector<uint8_t> pixels);
            
            /**
             * Initialize an image of size width * height with the given pixels.
             *
             * The array of pixel size must be equal to width * height.
             *
             * @param width Width if the image.
             * @param height Height of the image.
             * @param pixels Pixels composing the image.
             */
            Image8u1(uint32_t width, uint32_t height, uint8_t *pixels);
            
            
            ////////////////////////////////////////////////////////////////////
            ///////////////////////// DATA ACCESS //////////////////////////////
            ////////////////////////////////////////////////////////////////////
            
            /**
             * Return a reference to a pixel.
             *
             * @tparam T Integral type use to as index of the pixel.
             *
             * @param index Index of the pixel.
             *
             * @return A reference to the pixel corresponding to index.
             */
            template<typename T>
            [[nodiscard]] uint8_t &operator[](T index) {
                static_assert(std::is_integral<T>::value, "Integral type required.");
                return this->pixels[static_cast<uint32_t >(index)];
            }
            
            
            /**
             * Return the valu of a pixel.
             *
             * @tparam T Integral type use to as index of the pixel.
             *
             * @param index Index of the pixel.
             *
             * @return The value of the pixel corresponding to index.
             */
            template<typename T>
            [[nodiscard]] uint8_t at(T index) const {
                static_assert(std::is_integral<T>::value, "Integral type required.");
                return this->pixels.at(static_cast<uint32_t >(index));
            }
            
            
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
             * Save the image to path as a PNG.
             *
             * @param path Path to the output image.
             */
            void save(const char *path) const;
            
            /**
             * Create a copy of this image in dst, adding a border of size around the image.
             *
             * Border will be mirror reflection of the border elements (corresponding to BORDER_REFLECT on openCV).
             * Border's size cannot be greater the width or height.
             *
             * @param dst Destination of the bordered copy.
             * @param size Size of the border.
             */
            void borderedCopy(Image8u1 &dst, uint32_t size) const;
    };
}

#endif //BILATERALFILTER_IMAGE8U1_HPP
