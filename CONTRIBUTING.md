# The rules for contributions
In order to build an homogeneous code, here it is some rules.

## File naming
To avoid errors with `C` files, please use:
- `.hpp` for `C++` headers and
- `.cpp` for `C++` source files.

## Class naming, typenames, functions and variables
All stuff which belong to `Graphee` must be declared within,
the namespace `graphe::`.

Then specifically few rules:
- Class names: `MyRandomClass` start and distinguish words with **uppercase**, no `_` in class names.
- Typenames: `MyTypenameT` **same rules** than `class` and finishes with a `T`.
- Functions: `some_void_func` **only lowercase** and use `_` between words.
- Variables: **same rules** as `functions`.

## Braces and indentation
Some useful programs help us, here the `astyle` command to apply:
```
$ astyle --mode=c --style=ansi -s2 File1 [File2] [...]
```

## The future of C++\{11, 14\}
- I recommend to avoid old style `C` cast as `(const char*) my_array` and prefer 
`{static/dynamic/reinterpret}_cast<>`.
- A wonderful feature is `std::move`, do not hesitate to write move constructor
`T(T&&)`, and move assignment `T& operator=(T&&)`.
- Try to avoid low-level array allocation either with `malloc` xor `new`, use instead
`std::array` in order to minimize risks of memory leaks.

### References
- [The C++ Programming Language, *B. Stroustrup*, 4th Edition, Addison-Wesley, Upper Saddle River](http://www.stroustrup.com/4th.html)
- [Effective Modern C++, *S. Meyers*, 1st Edition, O'Reilly, Beijing](http://shop.oreilly.com/product/0636920033707.do)

**Thank you** :+1:

For any questions, comments, or collaborations, please use: **n.martin [at] qwantresearch [dot] com**
