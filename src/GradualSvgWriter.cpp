/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "GradualSvgWriter.hpp"
#include "ImmersedBoundaryMesh.hpp"

template <unsigned DIM>
GradualSvgWriter<DIM>::GradualSvgWriter()
        : AbstractCellBasedSimulationModifier<DIM>(),
          mSamplingMultiple(100u),
          mSvgSize(1600.0),
          mOutputDirectory(""),
          mSvgHeader(""),
          mSvgFooter("")
{
}

template <unsigned DIM>
GradualSvgWriter<DIM>::~GradualSvgWriter()
{
}

template <unsigned DIM>
void GradualSvgWriter<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    if (SimulationTime::Instance()->GetTimeStepsElapsed() % mSamplingMultiple == 0)
    {
        ImmersedBoundaryMesh<DIM, DIM>* p_mesh = static_cast<ImmersedBoundaryMesh<DIM, DIM>*>(&(rCellPopulation.rGetMesh()));

        // Get the number of time steps elapsed to use in the file name, and zero-pad it
        std::stringstream time;
        time << std::setfill('0') << std::setw(6) << SimulationTime::Instance()->GetTimeStepsElapsed();

        // Open the svg file for writing
        std::string full_file_name = "results_" + time.str() + ".svg";
        OutputFileHandler results_handler(mOutputDirectory, false);
        out_stream svg_file = results_handler.OpenOutputFile(full_file_name);

        (*svg_file) << mSvgHeader;

        // Add all nodes to the svg file
        double node_rad = p_mesh->GetAverageNodeSpacingOfElement(0, false) * 0.35 * mSvgSize;
        for (typename AbstractMesh<DIM, DIM>::NodeIterator it = p_mesh->GetNodeIteratorBegin();
             it != p_mesh->GetNodeIteratorEnd();
             ++it)
        {
            AddPointToSvgFile(svg_file, it->rGetLocation(), it->GetRegion(), node_rad);
        }

        // Add glyphs to the svg file
        unsigned num_elems_per_region = p_mesh->GetNumElements() / 3;
        double glyph_rad = 0.005 * mSvgSize;
        for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator it = p_mesh->GetElementIteratorBegin();
             it != p_mesh->GetElementIteratorEnd();
             ++it)
        {
            c_vector<double, 2> short_axis = p_mesh->GetShortAxisOfElement(it->GetIndex());
            short_axis = short_axis[0] < 0.0 ? short_axis : -short_axis;

            int angle = 90 + static_cast<int>(atan2(short_axis[0], short_axis[1]) * 180.0 / M_PI);

            AddGlyphToSvgFile(svg_file,
                              p_mesh->GetCentroidOfElement(it->GetIndex()),
                              static_cast<unsigned>(std::floor(it->GetIndex() / (double)num_elems_per_region)),
                              glyph_rad,
                              p_mesh->GetElongationShapeFactorOfElement(it->GetIndex()),
                              angle);
        }

        (*svg_file) << "<text x=\"" << 0.05 * mSvgSize << "\" y=\"" << (0.45 + 9.0 / 32.0) * mSvgSize << "\" "
                    << "font-size=\"30px\" font-family=\"monospace\">"
                    << "ts=" << SimulationTime::Instance()->GetTimeStepsElapsed()
                    << "</text>";

        (*svg_file) << mSvgFooter;

        svg_file->close();
    }
}

template <unsigned DIM>
void GradualSvgWriter<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;

    // Define colours
    std::string bg_col = "darkgray";
    std::string region0_col = "#990000"; // dark red
    std::string region1_col = "#cc0000"; // light red
    std::string region2_col = "#e68a00"; // dark orange
    std::string region3_col = "#ff9900"; // light orange
    std::string region4_col = "#006666"; // dark teal
    std::string region5_col = "#009999"; // light teal
    std::string region6_col = "#000099"; // dark blue
    std::string region7_col = "#0000cc"; // light blue
    std::string region8_col = "#FFFFFF"; // white
    std::string region9_col = "#000000"; // black
    std::stringstream header;

    // Set precision so as not to write out too many decimal places of uncompressed text
    header << std::setprecision(5);

    // Output svg header to file
    header << "<svg version=\"1.1\" baseProfile=\"full\" width=\""
           << mSvgSize << "px\" height=\"" << mSvgSize << "px\" "
           << "viewBox=\"0 0 " << mSvgSize << " " << mSvgSize << "\" "
           << "xmlns=\"http://www.w3.org/2000/svg\">"
           << "\n";

    // Add text/css style for elements
    header << "<style type=\"text/css\">"
           << "\n"
           << "\t.bg_rect{fill:" << bg_col << ";}"
           << "\n"
           << "\t.node_0{fill:" << region0_col << ";}"
           << "\n"
           << "\t.node_1{fill:" << region1_col << ";}"
           << "\n"
           << "\t.node_2{fill:" << region2_col << ";}"
           << "\n"
           << "\t.node_3{fill:" << region3_col << ";}"
           << "\n"
           << "\t.node_4{fill:" << region4_col << ";}"
           << "\n"
           << "\t.node_5{fill:" << region5_col << ";}"
           << "\n"
           << "\t.node_6{fill:" << region6_col << ";}"
           << "\n"
           << "\t.node_7{fill:" << region7_col << ";}"
           << "\n"
           << "\t.node_8{fill:" << region8_col << ";}"
           << "\n"
           << "\t.node_9{fill:" << region9_col << ";}"
           << "\n"
           << "</style>"
           << "\n";

    // Add background rectangle
    header << "<rect class=\"bg_rect\" width=\"" << mSvgSize << "\" height=\"" << mSvgSize << "\"/>"
           << "\n";

    mSvgHeader = header.str();
    mSvgFooter = "</svg>\n";
}

template <unsigned DIM>
void GradualSvgWriter<DIM>::AddPointToSvgFile(out_stream& rSvgFile, c_vector<double, DIM> location, unsigned region, double rad)
{
    double scaled_x = location[0] * mSvgSize;
    double scaled_y = (1.0 - location[1]) * mSvgSize;

    (*rSvgFile) << "<circle class=\"node_" << region << "\" "
                << "cx=\"" << scaled_x << "\" "
                << "cy=\"" << scaled_y << "\" "
                << "r=\"" << rad << "\"/>"
                << "\n";

    // Account for possible wrap-around of glyph in x
    if (scaled_x < rad)
    {
        (*rSvgFile) << "<circle class=\"node_" << region << "\" "
                    << "cx=\"" << scaled_x + mSvgSize << "\" "
                    << "cy=\"" << scaled_y << "\" "
                    << "r=\"" << rad << "\"/>"
                    << "\n";
    }
    else if (scaled_x > mSvgSize - rad)
    {
        (*rSvgFile) << "<circle class=\"node_" << region << "\" "
                    << "cx=\"" << scaled_x - mSvgSize << "\" "
                    << "cy=\"" << scaled_y << "\" "
                    << "r=\"" << rad << "\"/>"
                    << "\n";
    }

    // Account for possible wrap-around of glyph in y
    if (scaled_y < rad)
    {
        (*rSvgFile) << "<circle class=\"node_" << region << "\" "
                    << "cx=\"" << scaled_x << "\" "
                    << "cy=\"" << scaled_y + mSvgSize << "\" "
                    << "r=\"" << rad << "\"/>"
                    << "\n";
    }
    else if (scaled_y > mSvgSize - rad)
    {
        (*rSvgFile) << "<circle class=\"node_" << region << "\" "
                    << "cx=\"" << scaled_x << "\" "
                    << "cy=\"" << scaled_y - mSvgSize << "\" "
                    << "r=\"" << rad << "\"/>"
                    << "\n";
    }
}

template <unsigned DIM>
void GradualSvgWriter<DIM>::AddGlyphToSvgFile(out_stream& rSvgFile,
                                              c_vector<double, DIM> location,
                                              unsigned region,
                                              double rad,
                                              double elongation,
                                              int angle)
{
    double scaled_x = location[0] * mSvgSize;
    double scaled_y = (1.0 - location[1]) * mSvgSize;

    // Calculate the lum value (for glyph colour)
    unsigned col = 60 + static_cast<unsigned>(120.0 * fabs(0.5 - location[0]));

    double max_size = std::max(rad, rad * elongation);

    (*rSvgFile) << "<ellipse "
                << "transform=\"translate(" << scaled_x << " " << scaled_y
                << ") rotate(" << angle << ")\" "
                << "rx=\"" << rad << "\" "
                << "ry=\"" << elongation * rad << "\" "
                << "fill=\"rgb(" << col << "," << col << "," << col << ")\""
                << "/>"
                << "\n";

    // Account for possible wrap-around of glyph in x
    if (scaled_x < max_size)
    {
        (*rSvgFile) << "<ellipse "
                    << "transform=\"translate(" << scaled_x + mSvgSize << " " << scaled_y
                    << ") rotate(" << angle << ")\" "
                    << "rx=\"" << rad << "\" "
                    << "ry=\"" << elongation * rad << "\" "
                    << "fill=\"rgb(" << col << "," << col << "," << col << ")\""
                    << "/>"
                    << "\n";
    }
    else if (scaled_x > mSvgSize - max_size)
    {
        (*rSvgFile) << "<ellipse "
                    << "transform=\"translate(" << scaled_x - mSvgSize << " " << scaled_y
                    << ") rotate(" << angle << ")\" "
                    << "rx=\"" << rad << "\" "
                    << "ry=\"" << elongation * rad << "\" "
                    << "fill=\"rgb(" << col << "," << col << "," << col << ")\""
                    << "/>"
                    << "\n";
    }

    // Account for possible wrap-around of glyph in y
    if (scaled_y < max_size)
    {
        (*rSvgFile) << "<ellipse "
                    << "transform=\"translate(" << scaled_x << " " << scaled_y + mSvgSize
                    << ") rotate(" << angle << ")\" "
                    << "rx=\"" << rad << "\" "
                    << "ry=\"" << elongation * rad << "\" "
                    << "fill=\"rgb(" << col << "," << col << "," << col << ")\""
                    << "/>"
                    << "\n";
    }
    else if (scaled_y > mSvgSize - max_size)
    {
        (*rSvgFile) << "<ellipse "
                    << "transform=\"translate(" << scaled_x << " " << scaled_y - mSvgSize
                    << ") rotate(" << angle << ")\" "
                    << "rx=\"" << rad << "\" "
                    << "ry=\"" << elongation * rad << "\" "
                    << "fill=\"rgb(" << col << "," << col << "," << col << ")\""
                    << "/>"
                    << "\n";
    }
}

template <unsigned DIM>
void GradualSvgWriter<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template <unsigned DIM>
unsigned int GradualSvgWriter<DIM>::GetSamplingMultiple() const
{
    return mSamplingMultiple;
}

template <unsigned DIM>
void GradualSvgWriter<DIM>::SetSamplingMultiple(unsigned int samplingMultiple)
{
    mSamplingMultiple = samplingMultiple;
}

template <unsigned DIM>
double GradualSvgWriter<DIM>::GetSvgSize() const
{
    return mSvgSize;
}

template <unsigned DIM>
void GradualSvgWriter<DIM>::SetSvgSize(double svgSize)
{
    mSvgSize = svgSize;
}

// Explicit instantiation
template class GradualSvgWriter<1>;
template class GradualSvgWriter<2>;
template class GradualSvgWriter<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GradualSvgWriter)
