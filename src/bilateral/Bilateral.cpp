#include <cmath>
#include <cassert>

#include <bilateral/Bilateral.hpp>


namespace bilateral {
    
    static void gaussianKernel(std::vector<double> &kernel, double sigma, int32_t radius) {
        const double sigma2 = sigma * sigma;
        const double factor = 1. / (2. * M_PI * sigma2);
        const double divisor = 2. * sigma2;
        auto sGaussian = [factor, divisor](double value) {
            return factor * std::exp(-(value * value) / divisor);
        };
        
        const int32_t diameter = radius * 2 + 1;
        uint32_t index = 0;
        double sum = 0;
        
        kernel.resize(static_cast<uint32_t >(diameter * diameter));
        for (int32_t y = -radius; y <= radius; y++) {
            for (int32_t x = -radius; x <= radius; x++, index++) {
                kernel[index] = sGaussian(std::sqrt(y * y + x * x));
                sum += kernel[index];
            }
        }
        // Normalize the probability mass outside the kernel evenly to all pixels within the kernel
        for (double &x : kernel) {
            x /= sum;
        }
    }
    
    
    void Bilateral::naive(const Image8u1 &src, Image8u1 &dst, double sSpace, double sRange) {
        assert(sRange > 0 && "Color sigma must be greater than 0");
        assert(sSpace > 0 && "Space sigma must be greater than 0");
        
        const int32_t radius = static_cast<const int32_t>(std::round(sSpace * 1.5));
        const int32_t diameter = radius * 2 + 1;
        Image8u1 bordered;
        
        src.borderedCopy(bordered, static_cast<uint32_t>(radius));
        dst = Image8u1(src.getWidth(), src.getHeight());
        
        std::vector<double> sKernel;
        gaussianKernel(sKernel, sSpace, radius);
    
        const double sigma2 = sRange * sRange;
        const double factor = 1. / (2. * M_PI * sigma2);
        const double divisor = 2. * sigma2;
        auto rGaussian = [factor, divisor](double value) {
            return factor * std::exp(-(value * value) / divisor);
        };
        
        int32_t widthb = static_cast<int32_t>(bordered.getWidth());
        int32_t width = static_cast<int32_t>(dst.getWidth());
        int32_t height = static_cast<int32_t>(dst.getHeight());
        uint32_t indexk, center, index;
        double sum, wp, kmul;
        for (int32_t y = 0; y < height; y++) {
            for (int32_t x = 0; x < width; x++) {
                sum = 0;
                wp = 0;
                
                center = static_cast<uint32_t>((y + radius) * widthb + x + radius);
                for (int32_t ky = -radius; ky <= radius; ky++) {
                    for (int32_t kx = -radius; kx <= radius; kx++) {
                        indexk = static_cast<uint32_t>((ky + radius) * diameter + kx + radius);
                        index = static_cast<uint32_t>((y + ky + radius) * widthb + x + kx + radius);
                        kmul = sKernel[indexk] * rGaussian(std::abs(bordered[index] - bordered[center]));
                        wp += kmul;
                        sum += kmul * static_cast<double>(bordered[index]);
                    }
                }
                
                dst[y * width + x] = static_cast<uint8_t>(1 / wp * sum);
            }
        }
    }
}


