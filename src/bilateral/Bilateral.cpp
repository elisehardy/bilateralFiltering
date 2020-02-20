#include <bilateral/Bilateral.hpp>
#include <cmath>


namespace bilateral {
    
//    void Bilateral::naive(const Image8u1 &src, Image8u1 &dst, uint32_t diameter, float sColor, float sSpace) {
//        float colorCoef = -.5f / (sColor * sColor);
//        float spaceCoef = -.5f / (sSpace * sSpace);
//        uint32_t i, j, k, maxk, radius;
//
//        if (sColor <= 0) {
//            sColor = 1;
//        }
//        if (sSpace <= 0) {
//            sSpace = 1;
//        }
//
//        if (diameter <= 0) {
//            radius = static_cast<uint32_t>(std::round(sSpace * 1.5f));
//        }
//        else {
//            radius = diameter / 2;
//        }
//        radius = std::max(radius, 1u);
//        diameter = radius * 2 + 1;
//
//        Mat temp;
//        copyMakeBorder(src, temp, radius, radius, radius, radius, borderType);
//
//        vector<float> _color_weight(cn * 256);
//        vector<float> _space_weight(d *d);
//        vector<int> _space_ofs(d *d);
//        float *color_weight = &_color_weight[0];
//        float *space_weight = &_space_weight[0];
//        int *space_ofs = &_space_ofs[0];
//
//        // initialize color-related bilateral filter coefficients
//        for (i = 0; i < 256 * cn; i++) {
//            color_weight[i] = (float) std::exp(i * i * gauss_color_coeff);
//        }
//
//        // initialize space-related bilateral filter coefficients
//        for (i = -radius, maxk = 0; i <= radius; i++) {
//            for (j = -radius; j <= radius; j++) {
//                double r = std::sqrt((double) i * i + (double) j * j);
//                if (r > radius) {
//                    continue;
//                }
//                space_weight[maxk] = (float) std::exp(r * r * gauss_space_coeff);
//                space_ofs[maxk++] = (int) (i * temp.step + j * cn);
//            }
//        }
//
//        for (i = 0; i < size.height; i++) {
//            const uchar *sptr = temp.data + (i + radius) * temp.step + radius * cn;
//            uchar *dptr = dst.data + i * dst.step;
//
//            if (cn == 1) {
//                for (j = 0; j < size.width; j++) {
//                    float sum = 0, wsum = 0;
//                    int val0 = sptr[j];
//                    for (k = 0; k < maxk; k++) {
//                        int val = sptr[j + space_ofs[k]];
//                        float w = space_weight[k] * color_weight[std::abs(val - val0)];
//                        sum += val * w;
//                        wsum += w;
//                    }
//                    // overflow is not possible here => there is no need to use CV_CAST_8U
//                    dptr[j] = (uchar) cvRound(sum / wsum);
//                }
//            }
//            else {
//                assert(cn == 3);
//                for (j = 0; j < size.width * 3; j += 3) {
//                    float sum_b = 0, sum_g = 0, sum_r = 0, wsum = 0;
//                    int b0 = sptr[j], g0 = sptr[j + 1], r0 = sptr[j + 2];
//                    for (k = 0; k < maxk; k++) {
//                        const uchar *sptr_k = sptr + j + space_ofs[k];
//                        int b = sptr_k[0], g = sptr_k[1], r = sptr_k[2];
//                        float w = space_weight[k] * color_weight[std::abs(b - b0) +
//                                                                 std::abs(g - g0) + std::abs(r - r0)];
//                        sum_b += b * w;
//                        sum_g += g * w;
//                        sum_r += r * w;
//                        wsum += w;
//                    }
//                    wsum = 1.f / wsum;
//                    b0 = cvRound(sum_b * wsum);
//                    g0 = cvRound(sum_g * wsum);
//                    r0 = cvRound(sum_r * wsum);
//                    dptr[j] = (uchar) b0;
//                    dptr[j + 1] = (uchar) g0;
//                    dptr[j + 2] = (uchar) r0;
//                }
//            }
//        }
//    }
}



