# The input table: pastes from a spreadsheet, a French or an English CSV, a
# table copied from a PDF; units and decimal marks as people write them.

r = CLS.parse_table("Oxide\tmass %\nCaO\t63,7\nSiO2\t19,3\nTotal\t\n")
@test r.separator == "\t" && r.decimal_comma
@test [(x.formula, x.amount, x.unit) for x in r.rows] == [("CaO", 63.7, "%"), ("SiO2", 19.3, "%")]
@test r.skipped == ["Oxide\tmass %", "Total"]

r = CLS.parse_table("CaO;63,7;%\nNaCl;1,5;g\nK2SO4;0,2;mmol\n")
@test r.separator == ";"
@test [(x.formula, x.amount, x.unit) for x in r.rows] == [("CaO", 63.7, "%"), ("NaCl", 1.5, "g"), ("K2SO4", 0.2, "mmol")]

r = CLS.parse_table("CaO,63.7\nSiO2,19.3\n"; unit = "g")
@test r.separator == "," && !r.decimal_comma
@test [x.amount for x in r.rows] == [63.7, 19.3] && all(x -> x.unit == "g", r.rows)

# a table copied from a PDF: blanks between the cells, a decimal comma, a percent sign
r = CLS.parse_table("CaO   63,7 %\nAl2O3 5.7\nFe2O3 3,4%\n")
@test r.separator == "blanks"
@test [(x.formula, x.amount, x.unit) for x in r.rows] == [("CaO", 63.7, "%"), ("Al2O3", 5.7, "%"), ("Fe2O3", 3.4, "%")]

# thousands separated by a space, a non-breaking space; exponents; Windows line ends
r = CLS.parse_table("H2O;1 234,5;g\r\nNaCl;1 000;mg\r\nKCl;1e-3;mol/kg\r\n")
@test [x.amount for x in r.rows] == [1234.5, 1000.0, 1.0e-3]
@test [x.unit for x in r.rows] == ["g", "mg", "mol/kg"]

# units as people spell them
@test CLS.canonical_unit("wt.%") == "%"
@test CLS.canonical_unit(" mmol/l ") == "mmol/L"
@test CLS.canonical_unit("μmol") == "µmol"
@test CLS.canonical_unit("furlong") === nothing
r = CLS.parse_table("CaO\t1\tfurlong\n"; unit = "mol")
@test r.rows[1].unit == "mol"                 # an unknown unit cell leaves the default

@test CLS.parse_number("12,5"; decimal_comma = true) == 12.5
@test CLS.parse_number("12,5"; decimal_comma = false) === nothing
@test CLS.parse_number("abc"; decimal_comma = true) === nothing
@test isempty(CLS.parse_table("").rows)
@test_throws ArgumentError CLS.parse_table("CaO 1"; unit = "furlong")

# rows to and from JSON, and back to a spreadsheet
row = CLS.AmountRow("C3S", 60.0, "%")
@test CLS.amount_row(CLS.to_dict(row)) == row
@test CLS.plain(row) == CLS.to_dict(row)
@test CLS.table_text([row, CLS.AmountRow("C2S", 15, "g", "cement")]) == "C3S\t60.0\t%\nC2S\t15.0\tg"
@test_throws ArgumentError CLS.AmountRow("C3S", 1, "furlong")
@test_throws ArgumentError CLS.AmountRow("C3S", 1, "g", "greek")
@test_throws ArgumentError CLS.amount_row(3)
@test_throws ArgumentError CLS.amount_row(Dict("amount" => 1))
@test_throws ArgumentError CLS.amount_row(Dict("formula" => "x", "amount" => "1"))
@test_throws ArgumentError CLS.amount_row(Dict("formula" => "x", "amount" => 1, "unit" => 3))
@test_throws ArgumentError CLS.amount_row(Dict("formula" => "x", "amount" => 1, "notation" => 3))
