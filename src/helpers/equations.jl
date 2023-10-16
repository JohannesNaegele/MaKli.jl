macro equations(block)
    exprs = block.args

    # Filter for assignments
    assignments = filter(e -> isa(e, Expr) && e.head == :(=), exprs)

    # Extract variable names and their values as symbols with colons
    pairs = [Expr(:call, :(=>), :($(":$(a.args[1])")), a.args[2]) for a in assignments]

    return esc(Expr(:call, :Dict, pairs...))
end
