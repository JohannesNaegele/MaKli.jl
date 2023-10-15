macro equations(block)
    exprs = block.args

    # Filter for assignments
    assignments = filter(e -> isa(e, Expr) && e.head == :(=), exprs)

    # Extract variable names and their values
    pairs = [(a.args[1] => a.args[2]) for a in assignments]

    # Convert to a Dict expression
    return :(Dict($pairs...))
end