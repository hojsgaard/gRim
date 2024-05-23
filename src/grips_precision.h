#include <iostream>

// d = decimal places
template<int d> 
std::ostream& f(std::ostream &os){
    os.setf(std::ios_base::fixed, std::ios_base::floatfield); 
    os.precision(d); 
    return os; 
}

// w = width, d = decimal places
template<int w, int d> 
std::ostream& f(std::ostream &os){
    os.setf(std::ios_base::fixed, std::ios_base::floatfield); 
    os.precision(d); 
    os.width(w);
    return os; 
}

// d = decimal places
template<int d> 
std::ostream& e(std::ostream &os){
    os.setf(std::ios_base::scientific, std::ios_base::floatfield); 
    os.precision(d); 
    return os; 
}

// w = width, d = decimal places
template<int w, int d> 
std::ostream& e(std::ostream &os){
    os.setf(std::ios_base::scientific, std::ios_base::floatfield); 
    os.precision(d); 
    os.width(w);
    return os; 
}
