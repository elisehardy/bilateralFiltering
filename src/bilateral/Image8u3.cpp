#include <algorithm>

#include <bilateral/Image8u3.hpp>
#include <bilateral/lodepng.hpp>


namespace bilateral {
    
    Image8u3::Image8u3(const char *filename) {
        std::vector<uint8_t> raw;
    
        decodePNG(raw, this->width, this->height, filename);
        
        this->size = this->width * this->height;
        this->pixels.resize(this->size);
        // lodepng loads in RGBA, we need to convert it in RGB.
        for (uint32_t i = 0; i < this->size; i++) {
            this->pixels[i] = { raw[i * 4], raw[i * 4 + 1], raw[i * 4 + 2] };
        }
    }
    
    
    Image8u3::Image8u3(const std::string &filename) :
            Image8u3(filename.c_str()) {
    }
    
    
    Image8u3::Image8u3(uint32_t width, uint32_t height) :
            height(height), width(width) {
        this->size = height * width;
        this->pixels.resize(height * width);
    }
    
    
    Image8u3::Image8u3(uint32_t width, uint32_t height, std::vector<Pixel> pixels) :
            pixels(std::move(pixels)), height(height), width(width), size(height * width) {
    }
    
    
    Image8u3::Image8u3(uint32_t width, uint32_t height, Pixel *pixels) :
            pixels(pixels, pixels + width * height), height(height), width(width), size(height * width) {
    }
    
    
    uint32_t Image8u3::getWidth() const {
        return this->width;
    }
    
    
    uint32_t Image8u3::getHeight() const {
        return this->height;
    }
    
    
    uint32_t Image8u3::getSize() const {
        return this->size;
    }
    
    
    void Image8u3::save(const char *path) const {
        std::vector<uint8_t> raw(this->size * 4);
        
        for (uint32_t i = 0; i < this->size; i++) {
            raw[i * 4] = this->pixels[i].r;
            raw[i * 4 + 1] = this->pixels[i].g;
            raw[i * 4 + 2] = this->pixels[i].b;
            raw[i * 4 + 3] = 255;
        }
    
        encodePNG(path, raw, this->width, this->height);
    }
    
    
    void Image8u3::borderedCopy(Image8u3 &dst, uint32_t size) const {
        uint32_t width = this->width + size + size;
        uint32_t height = this->height + size + size;
        uint32_t dstY, srcY;
        
        dst = Image8u3(width, height);
        
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
}
