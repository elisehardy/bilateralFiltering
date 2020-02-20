

#include <bilateral/Image8u1.hpp>


using namespace bilateral;


int main(int argc, char **argv) {
    Image8u1 src = Image8u1(argv[1]), bordered;
    
    src.borderedCopy(bordered, 300);
    
    bordered.save("output.png");
    return 0;
}
