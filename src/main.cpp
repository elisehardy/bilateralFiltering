#include <getopt.h>

#include <bilateral/Bilateral.hpp>

#define DEFAULT_RANGE_SIGMA 64.

using namespace bilateral;


[[noreturn]] static void usageAndExit() {
    std::cout << std::endl;
    std::cout << "Usage:" << std::endl;
    std::cout << "\t  bilateralFilter [OPTIONS] -s SPACE_SIGMA -r RANGE_SIGMA INPUT_IMAGE" << std::endl << std::endl;
    
    std::cout << "Apply a bilateral filter with given sigma on a PNG image." << std::endl << std::endl;
    
    std::cout << "Mandatory arguments: " << std::endl << std::endl;
    
    std::cout << "\t -s SPACE_SIGMA" << std::endl;
    std::cout << "\t\t Sigma use in the gaussian function for the space weight," << std::endl;
    std::cout << "\t\t default to 2% of the input's diagonal" << std::endl << std::endl;
    
    std::cout << "\t -r RANGE_SIGMA" << std::endl;
    std::cout << "\t\t Sigma use in the gaussian function for the range weight," << std::endl;
    std::cout << "\t\t default to " << DEFAULT_RANGE_SIGMA << std::endl << std::endl;
    
    std::cout << "\t INPUT_IMAGE" << std::endl;
    std::cout << "\t\t Path to the image the bilateral filter will be applied on" << std::endl << std::endl;
    
    std::cout << "Optional arguments: " << std::endl << std::endl;
    
    std::cout << "\t -i ITERATION" << std::endl;
    std::cout << "\t\t Apply the bilateral filter ITERATION time (default to 1)" << std::endl << std::endl;
    
    std::cout << "\t -o" << std::endl;
    std::cout << "\t\t Path to the output image, default to 'output.png'" << std::endl << std::endl;
    
    std::cout << "\t -g" << std::endl;
    std::cout << "\t\t Read the image as a gray-scaled image, much faster than RGB." << std::endl << std::endl;
    
    std::cout << "\t -h" << std::endl;
    std::cout << "\t\t Displays this help" << std::endl << std::endl;
    exit(EXIT_SUCCESS);
}


int main(int argc, char **argv) {
    char c;
    bool gray = false;
    double sSigma = -1;
    double rSigma = -1;
    int32_t iteration = 1;
    std::string input;
    std::string output = "output.png";
    
    while ((c = static_cast<char>(getopt(argc, argv, "-s:-r:o:-i:gh"))) != -1) {
        switch (c) {
            case 's':
                try {
                    sSigma = std::stod(optarg);
                } catch (std::invalid_argument const &) {
                    std::cout << "Error: '" << optarg << "' is not a valid floating point number (see help with '-h')"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'r':
                try {
                    rSigma = std::stod(optarg);
                } catch (std::invalid_argument const &) {
                    std::cout << "Error: '" << optarg << "' is not a valid floating point number (see help with '-h')"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'i':
                try {
                    iteration = std::stoi(optarg);
                } catch (std::invalid_argument const &) {
                    std::cout << "Error: '" << optarg << "' is not a valid positive integer (see help with '-h')"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'o':
                output = optarg;
                break;
            case 'g':
                gray = true;
                break;
            case 'h':
                usageAndExit();
            case '?':
                if (optopt == 's' || optopt == 'r' || optopt == 'i' || optopt == 'o') {
                    std::cerr << "Option '-" << optopt << "' requires an argument (see help with '-h')" << std::endl;
                }
                else if (isprint(optopt)) {
                    std::cerr << "Unknown option '-" << optopt << "' (see help with '-h')" << std::endl;
                }
                else {
                    fprintf(stderr, "Unknown option character `\\x%x'.\n (see help with '-h')", optopt);
                }
                exit(EXIT_FAILURE);
            default:
                break;
        }
    }
    
    if (argc < 2 || optind > argc) {
        std::cerr << "Expected path to the image after options (see help with '-h')" << std::endl;
        exit(EXIT_FAILURE);
    }
    input = argv[optind - 1];
    
    if (gray) {
        Image8u1 src = Image8u1(input), dst;
        
        if (sSigma < 0) {
            uint32_t width = src.getWidth();
            uint32_t height = src.getHeight();
            sSigma = std::sqrt(width * width + height * height) * 0.02;
        }
        if (rSigma < 0) {
            rSigma = DEFAULT_RANGE_SIGMA;
        }
        
        for (int32_t i = 0; i < iteration; i++) {
            Bilateral::naive(src, dst, sSigma, rSigma);
            src = dst;
        }
        
        dst.save("output.png");
    }
    else {
        Image8u3 src = Image8u3(input), dst;
    
        if (sSigma < 0) {
            uint32_t width = src.getWidth();
            uint32_t height = src.getHeight();
            sSigma = std::sqrt(width * width + height * height) * 0.02;
        }
        if (rSigma < 0) {
            rSigma = DEFAULT_RANGE_SIGMA;
        }
        
        for (int32_t i = 0; i < iteration; i++) {
            Bilateral::naive(src, dst, sSigma, rSigma);
            src = dst;
        }
        
        dst.save("output.png");
    }
    
    return 0;
}
