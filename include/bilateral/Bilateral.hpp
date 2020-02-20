#ifndef BILATERALFILTER_BILATERAL_HPP
#define BILATERALFILTER_BILATERAL_HPP

#include <bilateral/Image8u1.hpp>


namespace bilateral {
    
    class Bilateral {
        
        public:
            
            Bilateral() = delete;
            
            static void naive(const Image8u1 &src, Image8u1 &dst, uint32_t diameter, float sColor, float sSpace);
    };
}

#endif // BILATERALFILTER_BILATERAL_HPP
