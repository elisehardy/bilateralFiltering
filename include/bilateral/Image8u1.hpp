#ifndef BILATERALFILTER_IMAGE8U1_HPP
#define BILATERALFILTER_IMAGE8U1_HPP

#include <string>
#include <utility>
#include <vector>
#include <iostream>


namespace bilateral {
    
    class Image8u1 {
        
        private:
            std::vector<uint8_t> pixels;
            uint32_t height = 0;
            uint32_t width = 0;
            uint32_t size = 0;
        
        public:
            
            ///////////////////////// CONSTRUCTORS /////////////////////////////////
            
            Image8u1() = default;
            
            
            explicit Image8u1(const char *filename);
            
            
            explicit Image8u1(const std::string &filename);
            
            
            Image8u1(uint32_t width, uint32_t height);
            
            
            Image8u1(uint32_t width, uint32_t height, std::vector<uint8_t> pixels);
            
            
            Image8u1(uint32_t width, uint32_t height, uint8_t *pixels);
            
            
            ///////////////////////// DATA ACCESS //////////////////////////////////
            
            template<typename T>
            [[nodiscard]] uint8_t &operator[](T index) {
                static_assert(std::is_integral<T>::value, "Integral type required.");
                return this->pixels[static_cast<uint32_t >(index)];
            }
            
            
            template<typename T>
            [[nodiscard]] uint8_t at(T index) const {
                static_assert(std::is_integral<T>::value, "Integral type required.");
                return this->pixels.at(static_cast<uint32_t >(index));
            }
            
            
            [[nodiscard]] uint32_t getWidth() const;
            
            
            [[nodiscard]] uint32_t getHeight() const;
            
            
            [[nodiscard]] uint32_t getSize() const;
            
            
            //////////////////////////// OTHERS ////////////////////////////////////
            
            void save(const char *filename) const;
            
            void borderedCopy(Image8u1 &dst, uint32_t size) const;
    };
}

#endif //BILATERALFILTER_IMAGE8U1_HPP
