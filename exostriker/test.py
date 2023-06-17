import cairosvg

def convert_svg_to_pdf(svg_file, pdf_file):
    try:
        cairosvg.svg2pdf(url=svg_file, write_to=pdf_file)
        print("Conversion completed successfully!")
    except Exception as e:
        print("Conversion failed:", str(e))

# Example usage
svg_file = "test.svg"
pdf_file = "output.pdf"

convert_svg_to_pdf(svg_file, pdf_file)

