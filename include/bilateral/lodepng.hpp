#ifndef BILATERALFILTER_LODEPNG_HPP
#define BILATERALFILTER_LODEPNG_HPP

#include <vector>
#include <cstdint>
#include <string>


namespace bilateral {
    
    void decodePNG(std::vector <uint8_t> &out, uint32_t &w, uint32_t &h, const std::string &filename);
    
    void encodePNG(const std::string &path, const std::vector <uint8_t> &in, uint32_t w, uint32_t h);
}

#endif //BILATERALFILTER_LODEPNG_HPP
