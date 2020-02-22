#ifndef BILATERALFILTER_BILATERAL_HPP
#define BILATERALFILTER_BILATERAL_HPP

#include <bilateral/Image8u1.hpp>


namespace bilateral {
    
    /**
     * Contains different implementation of the bilateral filter.
     */
    class Bilateral {
        
        public:
            
            Bilateral() = delete;
            
            /**
             * Compute in dst the result of a bilateral filter on src using different sigma for the space and range factor.
             *
             * @param src Image the bilateral filter will ba applied on.
             * @param dst Destination of the result of the bilateral filter on src.
             * @param sSpace Sigma use for the space weight.
             * @param sRange Sigma use for the range weight.
             */
            static void naive(const Image8u1 &src, Image8u1 &dst, double sSpace, double sRange);
    };
}

#endif // BILATERALFILTER_BILATERAL_HPP
