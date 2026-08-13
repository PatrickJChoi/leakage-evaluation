import os
import re
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, HRFlowable, PageBreak, KeepTogether
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER, TA_LEFT

def render_math(inner):
    # Unwrap \text{...} to its plain contents
    inner = re.sub(r"\\text\{([^{}]*)\}", r"\1", inner)
    inner = inner.replace("\\log", "log")
    inner = inner.replace("\\approx", "≈")
    inner = inner.replace("\\_", "_")
    # Drop any remaining LaTeX escape backslashes
    inner = inner.replace("\\", "")
    return inner

def clean_text(text):
    # Escape XML characters
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    # Replace bold and italic markdown markers with HTML tags
    text = re.compile(r"\*\*(.*?)\*\*").sub(r"<b>\1</b>", text)
    text = re.compile(r"\*(.*?)\*").sub(r"<i>\1</i>", text)
    # Clean LaTeX math delimiters: \( ... \), \[ ... \], and $ ... $
    text = re.compile(r"\\\((.*?)\\\)").sub(lambda m: f"<i>{render_math(m.group(1))}</i>", text)
    text = re.compile(r"\\\[(.*?)\\\]").sub(lambda m: f"<i>{render_math(m.group(1))}</i>", text)
    text = re.compile(r"\$([^\$]*)\$").sub(lambda m: f"<i>{render_math(m.group(1))}</i>", text)
    # Fallback: strip any stray $ or LaTeX escape characters that survived
    text = text.replace("$", "")
    text = re.sub(r"\\text\{([^{}]*)\}", r"\1", text)
    text = text.replace("\\", "")
    return text

def add_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont('Times-Roman', 9)
    canvas.setFillColor(colors.HexColor('#555555'))
    # Draw page number centered at the bottom of the page
    page_num = canvas.getPageNumber()
    canvas.drawCentredString(letter[0]/2.0, 40, f"Page {page_num}")
    canvas.restoreState()

def generate_pdf():
    _repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    preprint_path = os.path.join(_repo_root, "Leakage_Evaluation_Preprint.md")
    pdf_path = os.path.join(_repo_root, "Leakage_Evaluation_Preprint.pdf")
    
    with open(preprint_path, "r", encoding="utf-8") as f:
        text = f.read()

    # Split text into paragraphs/headers by double newlines
    blocks = [b.strip() for b in text.split("\n\n") if b.strip()]

    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=letter,
        leftMargin=79.2,   # 1.1 inches
        rightMargin=79.2,
        topMargin=79.2,
        bottomMargin=79.2
    )

    styles = getSampleStyleSheet()

    # Define academic document styles
    title_style = ParagraphStyle(
        'DocTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=16,
        leading=20,
        alignment=TA_CENTER,
        spaceAfter=15
    )
    
    author_style = ParagraphStyle(
        'DocAuthor',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=11,
        leading=14,
        alignment=TA_CENTER,
        spaceAfter=4
    )

    affiliation_style = ParagraphStyle(
        'DocAffiliation',
        parent=styles['Normal'],
        fontName='Helvetica-Oblique',
        fontSize=9.5,
        leading=12,
        alignment=TA_CENTER,
        spaceAfter=25
    )

    h1_style = ParagraphStyle(
        'SectionHeader',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=11,
        leading=14,
        alignment=TA_LEFT,
        spaceBefore=16,
        spaceAfter=6,
        keepWithNext=True
    )

    h2_style = ParagraphStyle(
        'SubSectionHeader',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=10.5,
        leading=13,
        alignment=TA_LEFT,
        spaceBefore=12,
        spaceAfter=4,
        keepWithNext=True
    )

    body_style = ParagraphStyle(
        'BodyJustified',
        parent=styles['Normal'],
        fontName='Times-Roman',
        fontSize=10,
        leading=14,
        alignment=TA_JUSTIFY,
        spaceAfter=8
    )

    abstract_body_style = ParagraphStyle(
        'AbstractBody',
        parent=styles['Normal'],
        fontName='Times-Roman',
        fontSize=9,
        leading=13,
        alignment=TA_JUSTIFY,
        leftIndent=24,
        rightIndent=24,
        spaceAfter=6
    )
    
    reference_style = ParagraphStyle(
        'ReferenceStyle',
        parent=styles['Normal'],
        fontName='Times-Roman',
        fontSize=8.5,
        leading=12,
        alignment=TA_LEFT,
        leftIndent=24,
        firstLineIndent=-24,
        spaceAfter=6
    )

    story = []
    
    # Process metadata blocks first
    title_text = clean_text(blocks[0].replace("# ", ""))
    story.append(Paragraph(title_text, title_style))
    
    author_text = clean_text(blocks[1])
    story.append(Paragraph(author_text, author_style))
    
    affil_text = clean_text(blocks[2])
    story.append(Paragraph(affil_text, affiliation_style))

    in_abstract = False
    in_references = False
    in_legends = False

    # Process remaining blocks
    for block in blocks[3:]:
        lines_in_block = block.split("\n")
        
        # Check if block starts with a header
        if block.startswith("# "):
            header_text = block.replace("# ", "").strip()
            
            if header_text.lower() == "abstract":
                in_abstract = True
                in_references = False
                in_legends = False
                story.append(Paragraph("<b>ABSTRACT</b>", ParagraphStyle('AbstractHeader', parent=title_style, fontSize=9.5, spaceAfter=8)))
            elif header_text.lower() == "references":
                in_abstract = False
                in_references = True
                in_legends = False
                story.append(Spacer(1, 10))
                story.append(Paragraph("REFERENCES", h1_style))
                story.append(HRFlowable(width="100%", thickness=0.75, color=colors.black, spaceBefore=2, spaceAfter=10))
            elif header_text.lower() == "figure legends":
                in_abstract = False
                in_references = False
                in_legends = True
                story.append(Spacer(1, 15))
                story.append(Paragraph("FIGURE LEGENDS", h1_style))
                story.append(HRFlowable(width="100%", thickness=0.75, color=colors.black, spaceBefore=2, spaceAfter=10))
            else:
                in_abstract = False
                in_references = False
                in_legends = False
                story.append(Spacer(1, 10))
                story.append(Paragraph(header_text.upper(), h1_style))
                story.append(HRFlowable(width="100%", thickness=0.75, color=colors.black, spaceBefore=2, spaceAfter=10))
                
        elif block.startswith("## "):
            header_text = block.replace("## ", "").strip()
            story.append(Paragraph(header_text, h2_style))
            
        elif block.startswith("### "):
            # Only the first line is the header; any remaining lines are body
            # text that happened to lack a blank-line separator in the source.
            header_text = lines_in_block[0].replace("### ", "").strip()
            story.append(Paragraph(f"<b><i>{header_text}</i></b>", ParagraphStyle('SubSubHeader', parent=h2_style, fontSize=10, spaceBefore=8)))
            body_lines = lines_in_block[1:]
            if body_lines:
                joined_block = " ".join([l.strip() for l in body_lines if l.strip()])
                cleaned = clean_text(joined_block)
                story.append(Paragraph(cleaned, body_style))

        else:
            # It's a standard text block
            if in_abstract:
                cleaned = clean_text(block)
                story.append(Paragraph(cleaned, abstract_body_style))
            elif in_references:
                # References can be multi-line in a single block
                for line in lines_in_block:
                    if line.strip():
                        cleaned = clean_text(line)
                        story.append(Paragraph(cleaned, reference_style))
            elif in_legends:
                # Figure legends
                for line in lines_in_block:
                    if line.strip():
                        cleaned = clean_text(line)
                        story.append(Paragraph(cleaned, ParagraphStyle('LegendText', parent=styles['Normal'], fontName='Times-Roman', fontSize=9, leading=13, spaceAfter=10, alignment=TA_JUSTIFY)))
            else:
                # Normal body text
                # Join lines if they were wrapped within a paragraph block
                joined_block = " ".join([l.strip() for l in lines_in_block])
                cleaned = clean_text(joined_block)
                story.append(Paragraph(cleaned, body_style))

    # Build the document
    doc.build(story, onFirstPage=add_footer, onLaterPages=add_footer)
    print(f"Preprint PDF generated successfully at: {pdf_path}")

if __name__ == "__main__":
    generate_pdf()
