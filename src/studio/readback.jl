# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# A script the studio wrote carries its model in a closing block comment, with
# the SHA-256 of the text above it. Reopening reads the model back exactly; a
# script edited by hand since (another fingerprint) is read too, but the studio
# does not take its path, so that saving never overwrites the hand edits.

const MODEL_TAG = "#= chemistrylab-studio-model v1"

# `#=` and `=#` open and close Julia block comments, which nest. Inside the
# JSON they can only occur within strings, where `#` (`#`) reads the same.
escape_block(json::AbstractString) = replace(json, "#=" => "\\u0023=", "=#" => "=\\u0023")

sha256_hex(text::AbstractString) = bytes2hex(SHA.sha256(text))

"""
    embed_model(body, m) -> String

`body` followed by the block that carries `m` and the fingerprint of `body`.
"""
function embed_model(body::AbstractString, m::StudioModel)
    json = escape_block(JSON.json(to_dict(m); pretty = true))
    return body * "\n" * MODEL_TAG * " sha256=" * sha256_hex(body) * "\n" * json * "\n=#\n"
end

"""
    ReadBack(model, edited, message)

What reading a script gives: the model (`nothing` when the file was not
written by the studio), whether its text was edited since, and a message for
the page.
"""
struct ReadBack
    model::Union{Nothing, StudioModel}
    edited::Bool
    message::String
end

"""
    read_script(source) -> ReadBack

Read the model a studio script carries.
"""
function read_script(source::AbstractString)
    k = findlast(MODEL_TAG, source)
    k === nothing && return ReadBack(nothing, false, "This file was not written by ChemistryLab Studio: it holds no model to read.")
    body = source[1:prevind(source, first(k))]
    endswith(body, '\n') && (body = chop(body))
    rest = source[nextind(source, last(k)):end]
    newline = findfirst('\n', rest)
    closing = findlast("=#", rest)
    (newline === nothing || closing === nothing) &&
        return ReadBack(nothing, false, "The model block at the end of this file is cut short.")
    m = match(r"sha256=([0-9a-f]{64})", rest[1:newline])
    json = rest[nextind(rest, newline):prevind(rest, first(closing))]
    model = model_from_dict(JSON.parse(json; dicttype = Dict{String, Any}))
    edited = m === nothing || m.captures[1] != sha256_hex(body)
    message = edited ?
        "This script was edited by hand since the studio wrote it: the studio shows its model, " *
        "and saving will ask for a new file name so that the edits are kept." : ""
    return ReadBack(model, edited, message)
end
