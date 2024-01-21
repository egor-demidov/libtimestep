//
// Created by egor on 1/19/24.
//

#ifndef INTEGRATORS_EXCEPTION_H
#define INTEGRATORS_EXCEPTION_H

#include <exception>
#include <string>

struct SizeMismatchException : std::exception {
    explicit SizeMismatchException(std::string const & source) :
            message("size mismatch in arguments passed to " + source) {}

    [[nodiscard]]
    const char * what() const noexcept override {
        return message.c_str();
    }

private:
    std::string message;
};

#endif //INTEGRATORS_EXCEPTION_H
