# python3 map_pdfs.py merged.pdf 2x1 ./../figs_QIIME2/qiime_vs_sampletype_family.pdf ./../figs/ampliconresequencing_vs_sampletype_family.pdf

from pypdf import PdfReader, PdfWriter, Transformation
import sys
import math

def merge_grid(pdf_paths, output_path, layout):
    cols, rows = map(int, layout.lower().split('x'))
    readers = [PdfReader(path) for path in pdf_paths]

    all_pages = []
    for reader in readers:
        all_pages.extend(reader.pages)

    writer = PdfWriter()

    pages_per_sheet = cols * rows
    total_pages = len(all_pages)
    sheets = math.ceil(total_pages / pages_per_sheet)

    for sheet_index in range(sheets):
        start = sheet_index * pages_per_sheet
        end = start + pages_per_sheet
        sheet_pages = all_pages[start:end]

        if not sheet_pages:
            continue

        widths = [p.mediabox.width for p in sheet_pages]
        heights = [p.mediabox.height for p in sheet_pages]

        max_width = max(widths)
        max_height = max(heights)

        new_width = cols * max_width
        new_height = rows * max_height

        # Create a blank new page with total size for the grid
        new_page = writer.add_blank_page(width=new_width, height=new_height)

        for idx, page in enumerate(sheet_pages):
            col_idx = idx % cols
            row_idx = idx // cols
            x_offset = col_idx * max_width
            y_offset = (rows - 1 - row_idx) * max_height

            # Merge the original page onto new_page with transformation
            new_page.merge_translated_page(page, tx=x_offset, ty=y_offset)

    with open(output_path, 'wb') as f_out:
        writer.write(f_out)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python merge_grid.py output.pdf NxM input1.pdf [input2.pdf ...]")
        print("Example: python merge_grid.py output.pdf 2x2 file1.pdf file2.pdf")
        sys.exit(1)

    output_file = sys.argv[1]
    layout = sys.argv[2]
    input_files = sys.argv[3:]

    merge_grid(input_files, output_file, layout)
