Rough guidelines to design principles:

- POD structs; no private members, no attached methods unless there's truly no better model.
Flow and lifetime/ownership issues are infinitely simpler this way.

- Prefer allocating and freeing blocks of memory,
rather than using complex smart ownership patterns for individual objects.
