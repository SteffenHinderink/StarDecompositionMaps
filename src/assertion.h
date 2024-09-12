#pragma once

#ifdef __FILE_NAME__
#define ASSERT_FILE_NAME __FILE_NAME__
#else
#define ASSERT_FILE_NAME __FILE__
#endif

#define ASSERT(expression)                                                                                                                                                            \
    do {                                                                                                                                                                              \
        if (!(expression)) {                                                                                                                                                          \
            std::cout << "Assertion (" << #expression << ") in file \"" << ASSERT_FILE_NAME << "\" function \"" << __FUNCTION__ << "\" line " << __LINE__ << " failed." << std::endl; \
            exit(0);                                                                                                                                                                  \
        }                                                                                                                                                                             \
    } while (false)
