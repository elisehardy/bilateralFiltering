#include <iostream>

#include <lodepng/lodepng.hpp>

#include <bilateral/lodepng.hpp>


namespace bilateral {
    
    void decodePNG(std::vector<uint8_t> &out, uint32_t &w, uint32_t &h, const std::string &filename) {
        uint32_t error = lodepng::decode(out, w, h, filename);
        if (error) {
            std::cerr << "Could not load image '" << filename << "' : " << lodepng_error_text(error) << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    
    void encodePNG(const std::string &path, const std::vector <uint8_t> &in, uint32_t w, uint32_t h) {
        uint32_t error = lodepng::encode(path, in, w, h);
        if (error) {
            std::cerr << "Could not save image to '" << path << "' : " << lodepng_error_text(error) << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}
