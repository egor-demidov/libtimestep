//
// Created by egor on 1/19/24.
//

#ifndef INTEGRATORS_EXCEPTION_H
#define INTEGRATORS_EXCEPTION_H

#include <exception>
#include <string>

// Exception thrown by system when the sizes of
// x0, v0, theta0, and omega0 buffers don't match
struct SizeMismatchException : std::exception {

    // Constructor takes one string - function/object that throws the exception
    explicit SizeMismatchException(std::string const & source) :
            message("size mismatch in arguments passed to " + source) {}

    // Returns the message with error description
    [[nodiscard]]
    const char * what() const noexcept override {
        return message.c_str();
    }

private:
    std::string message;
};

#endif //INTEGRATORS_EXCEPTION_H
