#include "phantom.h"

using namespace phantom;

/*
Helper function: Prints the parameters in a SEALContext.
*/
inline void print_parameters(const PhantomContext &context) {
  auto &context_data = context.get_context_data(0);
  /*
  Which scheme are we using?
  */
  std::string scheme_name;
  switch (context_data.parms().scheme()) {
    case phantom::scheme_type::bfv:
      scheme_name = "BFV";
      break;
    case phantom::scheme_type::ckks:
      scheme_name = "CKKS";
      break;
    case phantom::scheme_type::bgv:
      scheme_name = "BGV";
      break;
    default:
      throw std::invalid_argument("unsupported scheme");
  }
  std::cout << "/" << std::endl;
  std::cout << "| Library: PhantomFHE" << std::endl;
  std::cout << "| Encryption parameters:" << std::endl;
  std::cout << "|   scheme: " << scheme_name << std::endl;
  std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

  /*
  Print the size of the true (product) coefficient modulus.
  */
  std::cout << "|   coeff_modulus size: ";
  std::cout << context_data.total_coeff_modulus_bit_count() << " (";
  auto parms = context_data.parms();
  auto &coeff_modulus = parms.coeff_modulus();
  std::size_t coeff_modulus_size = coeff_modulus.size();
  for (std::size_t i = 0; i < coeff_modulus_size - 1; i++) {
    std::cout << coeff_modulus[i].bit_count() << " + ";
  }
  std::cout << coeff_modulus.back().bit_count();
  std::cout << ") bits" << std::endl;

  /*
  For the BFV scheme print the plain_modulus parameter.
  */
  if (context_data.parms().scheme() == phantom::scheme_type::bfv) {
    std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
  }

  std::cout << "\\" << std::endl;
}

/*
Helper function: Prints a vector of floating-point values.
*/
inline void print_vector(std::vector<cuDoubleComplex> vec, std::size_t print_size = 4, int prec = 3) {
  /*
  Save the formatting information for std::cout.
  */
  std::ios old_fmt(nullptr);
  old_fmt.copyfmt(std::cout);

  std::size_t slot_count = vec.size();

  std::cout << std::fixed << std::setprecision(prec);
  std::cout << std::endl;
  if (slot_count <= 2 * print_size) {
    std::cout << "    [";
    for (std::size_t i = 0; i < slot_count; i++) {
      std::cout << " " << vec[i].x << " + i * " << vec[i].y << ((i != slot_count - 1) ? "," : " ]\n");
    }
  } else {
    vec.resize(std::max(vec.size(), 2 * print_size));
    std::cout << "    [";
    for (std::size_t i = 0; i < print_size; i++) {
      std::cout << " " << vec[i].x << " + i * " << vec[i].y << ",";
    }
    if (vec.size() > 2 * print_size) {
      std::cout << " ...,";
    }
    for (std::size_t i = slot_count - print_size; i < slot_count; i++) {
      std::cout << " " << vec[i].x << " + i * " << vec[i].y << ((i != slot_count - 1) ? "," : " ]\n");
    }
  }
  std::cout << std::endl;

  /*
  Restore the old std::cout formatting.
  */
  std::cout.copyfmt(old_fmt);
}

inline bool near_equal(cuDoubleComplex x, cuDoubleComplex y, double tolerance) {
  if (std::abs(x.x - y.x) > tolerance) {
    return false;
  }

  return true;
}

inline bool near_equal(std::vector<double> x, std::vector<double> y, double tolerance) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  for (size_t i = 0; i < x.size(); i++) {
    if (std::abs(x[i] - y[i]) > tolerance) {
      return false;
    }
  }

  return true;
}
