#include <map>
#include <string>
#include <iostream>

int main(){
    std::map<int, std::string> m;
    m[1] = "uno";

    std::cout << m[1] << std::endl;

    return 0;
}