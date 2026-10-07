#!/usr/bin/env python3

"""
draw_ligand_table.py

Generate a publication-style table of ligand structures from a CSV file.

For each ligand, displays:
    - 2D chemical structure
    - Ligand name
    - SMILES string
    - Molecular formula
    - Molecular weight

Input CSV format:
    name,smiles
    PV-000123456789,CCOc1ccc(...)

A header is optional. If there is no header, the first column is assumed
to be the ligand name and the second column the SMILES string.

Example usage:
    python draw_ligand_table.py ligands.csv ligand_table.png

Optional:
    python draw_ligand_table.py ligands.csv ligand_table.png \
        --columns 5 \
        --cell-width 420 \
        --cell-height 520 \
        --dpi 300
"""

import argparse
import csv
import math
import os
from io import BytesIO

from PIL import Image, ImageDraw, ImageFont

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
from rdkit.Chem import rdDepictor


# ============================================================
# Font utilities
# ============================================================

def get_font(size, bold=False):
    """
    Try several commonly available TrueType fonts.

    DejaVu Sans is generally available on Linux systems and gives
    substantially better rendering than PIL's built-in bitmap font.
    """

    if bold:
        candidates = [
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
            "C:/Windows/Fonts/arialbd.ttf",
        ]
    else:
        candidates = [
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/dejavu/DejaVuSans.ttf",
            "C:/Windows/Fonts/arial.ttf",
        ]

    for path in candidates:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)

    # Pillow will use a default font if no TrueType font can be located.
    return ImageFont.load_default()


# ============================================================
# Text utilities
# ============================================================

def text_width(draw, text, font):
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0]


def text_height(draw, text, font):
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[3] - bbox[1]


def wrap_text_by_width(draw, text, font, max_width):
    """
    Wrap arbitrary text according to its rendered pixel width.

    This is particularly useful for SMILES because they normally contain
    no whitespace at which conventional text wrapping could occur.

    Newlines are inserted only for display; the original SMILES is not
    modified.
    """

    if not text:
        return [""]

    lines = []
    current = ""

    for char in text:
        candidate = current + char

        if text_width(draw, candidate, font) <= max_width:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = char

    if current:
        lines.append(current)

    return lines


def draw_centered_text(draw, text, y, font, cell_x, cell_width, fill="black"):
    width = text_width(draw, text, font)
    x = cell_x + (cell_width - width) / 2

    draw.text(
        (x, y),
        text,
        font=font,
        fill=fill,
    )


def draw_centered_multiline(
    draw,
    lines,
    start_y,
    font,
    cell_x,
    cell_width,
    line_spacing=3,
    fill="black",
):
    """
    Draw several individually centered lines.
    """

    y = start_y

    for line in lines:
        width = text_width(draw, line, font)
        x = cell_x + (cell_width - width) / 2

        draw.text(
            (x, y),
            line,
            font=font,
            fill=fill,
        )

        y += text_height(draw, line if line else "Ag", font) + line_spacing

    return y


# ============================================================
# CSV input
# ============================================================

def read_ligands(csv_file):
    """
    Read name/SMILES pairs from a CSV.

    Header is optional.

    Header names recognized include:
        name
        ligand
        ligand_name

    and:
        smiles
        smile
        canonical_smiles
        canonicalsmiles
    """

    with open(csv_file, "r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)

        rows = []

        for row in reader:
            # Skip empty rows
            if not row or all(not value.strip() for value in row):
                continue

            # Skip comments
            if row[0].strip().startswith("#"):
                continue

            rows.append(row)

    if not rows:
        raise ValueError("No ligand records found in CSV.")

    first = [x.strip().lower() for x in rows[0]]

    name_headers = {
        "name",
        "ligand",
        "ligand_name",
        "ligand name",
        "compound",
        "compound_name",
    }

    smiles_headers = {
        "smiles",
        "smile",
        "canonical_smiles",
        "canonical smiles",
        "canonicalsmiles",
    }

    has_header = (
        len(first) >= 2
        and first[0] in name_headers
        and first[1] in smiles_headers
    )

    if has_header:
        rows = rows[1:]

    ligands = []

    for row_number, row in enumerate(rows, start=2 if has_header else 1):

        if len(row) < 2:
            print(
                f"WARNING: row {row_number} contains fewer than two columns. "
                f"Skipping."
            )
            continue

        name = row[0].strip()
        smiles = row[1].strip()

        if not name or not smiles:
            print(
                f"WARNING: row {row_number} has an empty name or SMILES. "
                f"Skipping."
            )
            continue

        ligands.append((name, smiles))

    return ligands


# ============================================================
# Chemistry
# ============================================================

def prepare_molecule(smiles):
    """
    Parse a molecule and calculate properties.
    """

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    # Generate consistent 2D coordinates.
    rdDepictor.Compute2DCoords(mol)

    formula = rdMolDescriptors.CalcMolFormula(mol)

    # RDKit Descriptors.MolWt returns average molecular weight,
    # rather than exact monoisotopic mass.
    molecular_weight = Descriptors.MolWt(mol)

    return mol, formula, molecular_weight


def render_molecule(mol, width, height):
    """
    Render an RDKit molecule into a Pillow image.
    """

    # Use MolToImage for a clean white-background depiction.
    image = Draw.MolToImage(
        mol,
        size=(width, height),
        kekulize=True,
        wedgeBonds=True,
        fitImage=True,
    )

    if image.mode != "RGB":
        image = image.convert("RGB")

    return image


# ============================================================
# Main table drawing
# ============================================================

def make_ligand_table(
    ligands,
    output_file,
    columns=5,
    cell_width=420,
    cell_height=520,
    structure_height=300,
    margin=20,
    name_font_size=22,
    smiles_font_size=16,
    property_font_size=18,
    grid_width=2,
    dpi=300,
):
    """
    Generate the ligand table.
    """

    if not ligands:
        raise ValueError("No ligands to draw.")

    rows = math.ceil(len(ligands) / columns)

    image_width = columns * cell_width + 2 * margin
    image_height = rows * cell_height + 2 * margin

    canvas = Image.new(
        "RGB",
        (image_width, image_height),
        "white",
    )

    draw = ImageDraw.Draw(canvas)

    name_font = get_font(name_font_size, bold=True)
    smiles_font = get_font(smiles_font_size, bold=False)
    property_font = get_font(property_font_size, bold=False)

    # ---------------------------------------------------------
    # Draw every table cell
    # ---------------------------------------------------------

    for index in range(rows * columns):

        row = index // columns
        col = index % columns

        x0 = margin + col * cell_width
        y0 = margin + row * cell_height

        x1 = x0 + cell_width
        y1 = y0 + cell_height

        # Cell border
        draw.rectangle(
            [x0, y0, x1, y1],
            outline="black",
            width=grid_width,
        )

        # Empty table cells at the end are intentionally left blank.
        if index >= len(ligands):
            continue

        ligand_name, smiles = ligands[index]

        result = prepare_molecule(smiles)

        # -----------------------------------------------------
        # Handle invalid SMILES
        # -----------------------------------------------------

        if result is None:

            draw_centered_text(
                draw,
                ligand_name,
                y0 + 30,
                name_font,
                x0,
                cell_width,
            )

            draw_centered_text(
                draw,
                "INVALID SMILES",
                y0 + 90,
                property_font,
                x0,
                cell_width,
                fill="red",
            )

            smiles_lines = wrap_text_by_width(
                draw,
                smiles,
                smiles_font,
                cell_width - 30,
            )

            draw_centered_multiline(
                draw,
                smiles_lines,
                y0 + 140,
                smiles_font,
                x0,
                cell_width,
            )

            continue

        mol, formula, molecular_weight = result

        # -----------------------------------------------------
        # Molecule depiction
        # -----------------------------------------------------

        molecule_margin_x = 15
        molecule_margin_y = 10

        mol_width = cell_width - 2 * molecule_margin_x
        mol_height = structure_height - 2 * molecule_margin_y

        mol_image = render_molecule(
            mol,
            mol_width,
            mol_height,
        )

        paste_x = x0 + molecule_margin_x
        paste_y = y0 + molecule_margin_y

        canvas.paste(
            mol_image,
            (paste_x, paste_y),
        )

        # -----------------------------------------------------
        # Ligand name
        # -----------------------------------------------------

        text_y = y0 + structure_height + 3

        # Names can wrap too if exceptionally long.
        name_lines = wrap_text_by_width(
            draw,
            ligand_name,
            name_font,
            cell_width - 24,
        )

        text_y = draw_centered_multiline(
            draw,
            name_lines,
            text_y,
            name_font,
            x0,
            cell_width,
            line_spacing=2,
        )

        text_y += 5

        # -----------------------------------------------------
        # SMILES
        # -----------------------------------------------------

        smiles_lines = wrap_text_by_width(
            draw,
            smiles,
            smiles_font,
            cell_width - 24,
        )

        text_y = draw_centered_multiline(
            draw,
            smiles_lines,
            text_y,
            smiles_font,
            x0,
            cell_width,
            line_spacing=1,
        )

        text_y += 7

        # -----------------------------------------------------
        # Formula
        # -----------------------------------------------------

        draw_centered_text(
            draw,
            formula,
            text_y,
            property_font,
            x0,
            cell_width,
        )

        text_y += text_height(draw, formula, property_font) + 5

        # -----------------------------------------------------
        # Molecular weight
        # -----------------------------------------------------

        mw_string = f"{molecular_weight:.2f} g/mol"

        draw_centered_text(
            draw,
            mw_string,
            text_y,
            property_font,
            x0,
            cell_width,
        )

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------

    extension = os.path.splitext(output_file)[1].lower()

    if extension in [".jpg", ".jpeg"]:
        canvas.save(
            output_file,
            quality=95,
            dpi=(dpi, dpi),
        )
    else:
        canvas.save(
            output_file,
            dpi=(dpi, dpi),
        )

    print()
    print(f"Created ligand table:")
    print(f"  {output_file}")
    print()
    print(f"Ligands: {len(ligands)}")
    print(f"Rows:     {rows}")
    print(f"Columns:  {columns}")
    print(f"Image:    {image_width} x {image_height} px")
    print(f"DPI:      {dpi}")


# ============================================================
# Command line interface
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate a table containing 2D ligand structures, names, "
            "SMILES strings, molecular formulas, and molecular weights."
        )
    )

    parser.add_argument(
        "input_csv",
        help="CSV containing ligand name followed by SMILES",
    )

    parser.add_argument(
        "output_image",
        help="Output image, e.g. ligand_table.png",
    )

    parser.add_argument(
        "--columns",
        type=int,
        default=5,
        help="Maximum number of ligands per row (default: 5)",
    )

    parser.add_argument(
        "--cell-width",
        type=int,
        default=420,
        help="Width of each table cell in pixels (default: 420)",
    )

    parser.add_argument(
        "--cell-height",
        type=int,
        default=520,
        help="Height of each table cell in pixels (default: 520)",
    )

    parser.add_argument(
        "--structure-height",
        type=int,
        default=300,
        help="Vertical space allocated to molecule structure (default: 300)",
    )

    parser.add_argument(
        "--margin",
        type=int,
        default=20,
        help="Outer image margin in pixels (default: 20)",
    )

    parser.add_argument(
        "--name-font-size",
        type=int,
        default=22,
        help="Ligand-name font size (default: 22)",
    )

    parser.add_argument(
        "--smiles-font-size",
        type=int,
        default=16,
        help="SMILES font size (default: 16)",
    )

    parser.add_argument(
        "--property-font-size",
        type=int,
        default=18,
        help="Formula/MW font size (default: 18)",
    )

    parser.add_argument(
        "--grid-width",
        type=int,
        default=2,
        help="Table border width (default: 2)",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output image DPI metadata (default: 300)",
    )

    args = parser.parse_args()

    ligands = read_ligands(args.input_csv)

    print(f"Read {len(ligands)} ligand(s) from {args.input_csv}")

    make_ligand_table(
        ligands=ligands,
        output_file=args.output_image,
        columns=args.columns,
        cell_width=args.cell_width,
        cell_height=args.cell_height,
        structure_height=args.structure_height,
        margin=args.margin,
        name_font_size=args.name_font_size,
        smiles_font_size=args.smiles_font_size,
        property_font_size=args.property_font_size,
        grid_width=args.grid_width,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
