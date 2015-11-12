###############################################################################
# Check for the presence of std::exception_ptr and std::rethrow
macro(NMZ_CHECK_FOR_EXCEPTION_PTR)

    # We need to check whether the compiler supports the rethrowing mechanism
    include(CheckCXXSourceRuns)
    # at least the c++11 flag needs to be set
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
    check_cxx_source_runs("
        #include <iostream>
        #include <exception>
        #include <stdexcept>

        int main () {
            std::exception_ptr tmp_exception;
            try {
                throw std::overflow_error(\"some overflow exception\");
            } catch(const std::exception& e) {
                tmp_exception = std::current_exception();
                std::cout << \"exception caught, but continuing...\" << std::endl;
            }

            std::cout << \"(after exception)\" << std::endl;

            try {
                if (tmp_exception != std::exception_ptr()) std::rethrow_exception(tmp_exception);
            } catch (const std::exception& e) {
                std::cout << \"exception caught again \" << e.what() << std::endl;
            }
            return 0;
        }
"
        HAVE_EXCEPTION_PTR)

endmacro(NMZ_CHECK_FOR_EXCEPTION_PTR)
