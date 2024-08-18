### C++ Style Guide

This style guide outlines the coding standards and best practices for developing and maintaining a C++ codebase. It is intended to ensure consistency, readability, and maintainability across the project.

---

## Table of Contents
1. [File Structure](#file-structure)
2. [Naming Conventions](#naming-conventions)
3. [Commenting and Documentation](#commenting-and-documentation)
4. [Code Formatting](#code-formatting)
5. [Classes and Functions](#classes-and-functions)
6. [Memory Management](#memory-management)
7. [Error Handling](#error-handling)
8. [Best Practices](#best-practices)

---

## File Structure

- **File Names**: Use `PascalCase` for file names (e.g., `MyClass.cpp`, `Utils.h`).
- **Header Files**: Place class declarations and function prototypes in `.h` files.
- **Source Files**: Place class definitions and function implementations in `.cpp` files.
- **Include Guards**: Use traditional include guards (`#ifndef`, `#define`, `#endif`) in header files.

```cpp
#ifndef MY_CLASS_H
#define MY_CLASS_H

// Code

#endif // MY_CLASS_H
```

## Naming Conventions

- **Variables**: Use `camelCase` for variables (e.g., `totalCount`, `userName`).
- **Member variables**: Use `camelCase` followed by an underscore for member variables (e.g., `totalCount_`, `userName_`). N.B.: for structs, do not add the underscore.
- **Functions**: Use `camelCase` for function names (e.g., `getUserName()`, `calculateSum()`).
- **Classes**: Use `PascalCase` for class names (e.g., `UserAccount`, `DataManager`).
- **Constants**: Use `ALL_CAPS_WITH_UNDERSCORES` for constants (e.g., `MAX_BUFFER_SIZE`).
- **Namespaces**: Use `lowercase` for namespaces (e.g., `namespace utils { ... }`).

## Commenting and Documentation

- **Block Comments**: Use block comments (`/* ... */`) for large, detailed descriptions.
- **Inline Comments**: Use inline comments (`//`) for brief explanations or notes.
- **Javadoc Style Comments**: Use Javadoc-style comments for function and class documentation.

### Function Documentation Example

```cpp
/**
 * Calculates the factorial of a number.
 *
 * @param n The number to calculate the factorial for.
 * @return The factorial of the given number.
 */
int factorial(int n);
```

### Class Documentation Example

```cpp
/**
 * Represents a user account in the system.
 */
class UserAccount {
public:
    /**
     * Constructor that initializes the user account with a name.
     *
     * @param name The name of the user.
     */
    UserAccount(const std::string& name);

    /**
     * Gets the user's name.
     *
     * @return The name of the user.
     */
    std::string getName() const;

private:
    std::string name_;
};
```

## Code Formatting

We use the [Chromium C++](https://chromium.googlesource.com/chromium/src/+/refs/heads/main/styleguide/c++/c++.md) formatting. Note that most files have not been formatted in this manner yet. This reformatting will be applied step-wise.

Below are the most important rules:

- **Indentation**: Use 2 spaces for indentation. Avoid tabs.
- **Line Length**: Limit lines to 80 characters where possible.
- **Braces**: 
  - Place opening braces on the same line as the control statement or function declaration.
  - Place closing braces on a new line.

```cpp
if (condition) {
    // Code block
} else {
    // Code block
}
```

- **Spaces**:
  - Use spaces around operators (e.g., `a + b`, `x = y`).
  - Do not add spaces after `(`, `[` or before `)`, `]`.

```cpp
int result = (a + b) * (c - d);
```

- **Blank Lines**: 
  - Use blank lines to separate different sections of code (e.g., between functions, after variable declarations).
  - Avoid excessive blank lines.

## Classes and Functions

- **Class Structure**:
  - Order class members as: public, protected, private.
  - Group related member functions together.

- **Function Length**:
  - Keep functions short and focused. Ideally, functions should not exceed 40 lines of code.
  - If a function becomes too long, consider breaking it into smaller, more manageable functions.

- **Parameter Passing**:
  - Pass objects by reference or pointer when possible to avoid unnecessary copying.
  - Use `const` where appropriate to indicate that a parameter should not be modified.

```cpp
void process(const std::vector<int>& data);
```

## Memory Management

- **Smart Pointers**: Prefer using smart pointers (`std::unique_ptr`, `std::shared_ptr`) over raw pointers to manage dynamic memory.
- **RAII**: Follow the RAII (Resource Acquisition Is Initialization) principle to ensure that resources are properly cleaned up.
- **Avoid Manual Memory Management**: Avoid using `new` and `delete` directly. Use containers and smart pointers instead.

## Error Handling

- **Exceptions**: Use exceptions for error handling instead of return codes. Throw exceptions when an error occurs that cannot be handled locally.
- **Error Messages**: Provide clear and informative error messages when throwing exceptions.
- **Assertions**: Use `assert()` for conditions that should never occur. Assertions should be used to catch programming errors, not runtime errors.

```cpp
void processData(int value) {
    assert(value >= 0); // Ensure the value is non-negative
    // Process the data
}
```

## Best Practices

- **Code Review**: All code should be reviewed by at least one other team member before being merged into the main branch.
- **Unit Testing**: Write unit tests for all critical functions. Aim for high test coverage.
- **Consistent Style**: Follow this style guide consistently across the entire codebase. Consistency is more important than individual preferences.
- **Avoid Global Variables**: Minimize the use of global variables. Use singletons or dependency injection when global state is needed.
- **Encapsulation**: Keep class members private unless they need to be accessed externally. Provide public getter and setter functions as needed.
- **Use STL**: Prefer using the Standard Template Library (STL) over custom data structures or algorithms where applicable.

---

This guide serves as a reference for writing clear, maintainable, and consistent C++ code. Adherence to these guidelines will improve collaboration and code quality across the project.