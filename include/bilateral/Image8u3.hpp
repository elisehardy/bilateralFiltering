#ifndef BILATERALFILTER_IMAGE8U3_HPP
#define BILATERALFILTER_IMAGE8U3_HPP

#include <string>
#include <utility>
#include <vector>
#include <iostream>

#include <glm/glm.hpp>


namespace bilateral {
    
    typedef glm::vec<3, uint8_t> Pixel;
    
    
    
    class Image8u3 {
        
        private:
            std::vector<Pixel> pixels;
            uint32_t height = 0;
            uint32_t width = 0;
            uint32_t size = 0;
        
        public:
            
            ///////////////////////// CONSTRUCTORS /////////////////////////////////
            
            Image8u3() = default;
            
            explicit Image8u3(const char *filename);
            
            explicit Image8u3(const std::string &filename);
            
            Image8u3(uint32_t width, uint32_t height);
            
            Image8u3(uint32_t width, uint32_t height, std::vector<Pixel> pixels);
            
            Image8u3(uint32_t width, uint32_t height, Pixel *pixels);
            
            
            ///////////////////////// DATA ACCESS //////////////////////////////////
            
            template<typename T>
            [[nodiscard]] Pixel &operator[](T index) {
                static_assert(std::is_integral<T>::value, "Integral type required.");
                return this->pixels[static_cast<uint32_t >(index)];
            }
            
            
            template<typename T>
            [[nodiscard]] Pixel at(T index) const {
                static_assert(std::is_integral<T>::value, "Integral type required.");
                return this->pixels.at(static_cast<uint32_t >(index));
            }
            
            
            [[nodiscard]] uint32_t getWidth() const;
            
            [[nodiscard]] uint32_t getHeight() const;
            
            [[nodiscard]] uint32_t getSize() const;
            
            
            //////////////////////////// OTHERS ////////////////////////////////////
            
            void save(const char *path) const;
            
            void borderedCopy(Image8u3 &dst, uint32_t size) const;
    };
}

#endif //BILATERALFILTER_IMAGE8U3_HPP
