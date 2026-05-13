Rough guidelines to design principles:

- Prefer simple aggregate structs; no private members, user-defined constructors,
or attached methods unless there's a strong reason.
Lifetime/ownership, and use in bindings, are simplified for doing so.

- Prefer allocating and freeing blocks of memory,
rather than using complex smart ownership patterns for individual objects.
