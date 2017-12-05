/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "ThreeRegionSvgWriter.hpp"
#include "ImmersedBoundaryMesh.hpp"

template <unsigned DIM>
ThreeRegionSvgWriter<DIM>::ThreeRegionSvgWriter()
        : AbstractCellBasedSimulationModifier<DIM>(),
          mSamplingMultiple(100u),
          mSvgSize(1600.0),
          mOutputDirectory(""),
          mSvgHeader(""),
          mSvgFooter(""),
          mRegionSizes()
{
}

template <unsigned DIM>
void ThreeRegionSvgWriter<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    if (SimulationTime::Instance()->GetTimeStepsElapsed() % mSamplingMultiple == 0)
    {
        auto p_mesh = static_cast<ImmersedBoundaryMesh<DIM, DIM> *>(&(rCellPopulation.rGetMesh()));

        // Get the number of time steps elapsed to use in the file name, and zero-pad it
        std::stringstream time;
        time << std::setfill('0') << std::setw(6) << SimulationTime::Instance()->GetTimeStepsElapsed();

        // Open the svg file for writing
        std::string full_file_name = "results_" + time.str() + ".svg";
        OutputFileHandler results_handler(mOutputDirectory, false);
        out_stream svg_file = results_handler.OpenOutputFile(full_file_name);

        (*svg_file) << mSvgHeader;

        // Add all nodes to the svg file
        const double node_rad = p_mesh->GetAverageNodeSpacingOfElement(0, false) * 0.35 * mSvgSize;
        for (typename AbstractMesh<DIM, DIM>::NodeIterator it = p_mesh->GetNodeIteratorBegin();
             it != p_mesh->GetNodeIteratorEnd();
             ++it)
        {
            AddPointToSvgFile(svg_file, it->rGetLocation(), it->GetRegion(), node_rad);
        }

        // Add glyphs to the svg file
        const double glyph_rad = 0.005 * mSvgSize;
        for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator it = p_mesh->GetElementIteratorBegin();
             it != p_mesh->GetElementIteratorEnd();
             ++it)
        {
            c_vector<double, 2> short_axis = p_mesh->GetShortAxisOfElement(it->GetIndex());
            short_axis = short_axis[0] < 0.0 ? short_axis : -short_axis;

            unsigned region = 2u;
            if (it->GetIndex() < mRegionSizes[0])
            {
                region = 0u;
            }
            else if (it->GetIndex() < mRegionSizes[0] + mRegionSizes[1])
            {
                region = 1u;
            }

            int angle = 90 + static_cast<int>(atan2(short_axis[0], short_axis[1]) * 180.0 / M_PI);

            AddGlyphToSvgFile(svg_file,
                              p_mesh->GetCentroidOfElement(it->GetIndex()),
                              region,
                              glyph_rad,
                              p_mesh->GetElongationShapeFactorOfElement(it->GetIndex()),
                              angle);

            // Add large black points to visualise corners
            for (const auto& corner : it->rGetCornerNodes())
            {
                AddPointToSvgFile(svg_file, corner->rGetLocation(), 8u, 1.5 * node_rad);
            }
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
void ThreeRegionSvgWriter<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;

    // Define colours
    const std::string bg_col = "darkgray";
    const std::string region0_col = "#990000"; // dark red
    const std::string region1_col = "#cc0000"; // light red
    const std::string region2_col = "#e68a00"; // dark orange
    const std::string region3_col = "#ff9900"; // light orange
    const std::string region4_col = "#006666"; // dark teal
    const std::string region5_col = "#009999"; // light teal
    const std::string region6_col = "#000099"; // dark blue
    const std::string region7_col = "#0000cc"; // light blue
    const std::string region8_col = "#FFFFFF"; // white
    const std::string region9_col = "#000000"; // black
    const std::string glyph0_col = "DarkRed"; // white
    const std::string glyph1_col = "DarkBlue"; // white
    const std::string glyph2_col = "DarkGreen"; // white

    std::stringstream header;

    // Set precision so as not to write out too many decimal places of uncompressed text
    header << std::setprecision(5);

    //
    double aspect_ratio = 16.0 / 9.0;
    double height = mSvgSize / aspect_ratio;
    double offset = 0.5 * (mSvgSize - height);

    // Output svg header to file
    header << "<svg version=\"1.1\" baseProfile=\"full\" width=\""
           << mSvgSize <<"px\" height=\"" << height << "px\" "
           << "viewBox=\"0 " << offset << " " << mSvgSize << " " << height << "\" "
           << "xmlns=\"http://www.w3.org/2000/svg\">" << "\n";

    // Add text/css style for elements
    header << "<style type=\"text/css\">" << "\n"
           << "\t.bg_rect{fill:" << bg_col << ";}" << "\n"
           << "\t.node_0{fill:" << region0_col << ";}" << "\n"
           << "\t.node_1{fill:" << region1_col << ";}" << "\n"
           << "\t.node_2{fill:" << region2_col << ";}" << "\n"
           << "\t.node_3{fill:" << region3_col << ";}" << "\n"
           << "\t.node_4{fill:" << region4_col << ";}" << "\n"
           << "\t.node_5{fill:" << region5_col << ";}" << "\n"
           << "\t.node_6{fill:" << region6_col << ";}" << "\n"
           << "\t.node_7{fill:" << region7_col << ";}" << "\n"
           << "\t.node_8{fill:" << region8_col << ";}" << "\n"
           << "\t.node_9{fill:" << region9_col << ";}" << "\n"
           << "\t.glyph_0{fill:" << glyph0_col << ";}" << "\n"
           << "\t.glyph_1{fill:" << glyph1_col << ";}" << "\n"
           << "\t.glyph_2{fill:" << glyph2_col << ";}" << "\n"
           << "</style>" << "\n";

    // Add background rectangle
    header << "<rect class=\"bg_rect\" width=\"" << mSvgSize << "\" height=\"" << mSvgSize << "\"/>" << "\n";

    mSvgHeader = header.str();
    mSvgFooter = "</svg>\n";

    ImmersedBoundaryMesh<DIM, DIM>* p_mesh = static_cast<ImmersedBoundaryMesh<DIM, DIM> *>(&(rCellPopulation.rGetMesh()));
    if (std::accumulate(mRegionSizes.begin(), mRegionSizes.end(), 0u) != p_mesh->GetNumElements())
    {
        EXCEPTION("mRegionSizes must be set to the correct number of elements");
    }
}

template <unsigned DIM>
void ThreeRegionSvgWriter<DIM>::AddPointToSvgFile(out_stream& rSvgFile, c_vector<double, DIM> location, unsigned region,
                                                  double rad) const noexcept
{
    const double scaled_x = location[0] * mSvgSize;
    const double scaled_y = (1.0 - location[1]) * mSvgSize;

    (*rSvgFile) << "<circle class=\"node_" << region << "\" "
                << "cx=\"" << scaled_x << "\" "
                << "cy=\"" << scaled_y << "\" "
                << "r=\"" << rad << "\"/>" << "\n";

    // Account for possible wrap-around of glyph in x
    if (scaled_x < rad)
    {
        (*rSvgFile) << "<circle class=\"node_" << region << "\" "
                    << "cx=\"" << scaled_x + mSvgSize << "\" "
                    << "cy=\"" << scaled_y << "\" "
                    << "r=\"" << rad << "\"/>" << "\n";
    }
    else if (scaled_x > mSvgSize - rad)
    {
        (*rSvgFile) << "<circle class=\"node_" << region << "\" "
                    << "cx=\"" << scaled_x - mSvgSize << "\" "
                    << "cy=\"" << scaled_y << "\" "
                    << "r=\"" << rad << "\"/>" << "\n";
    }

    // Account for possible wrap-around of glyph in y
    if (scaled_y < rad)
    {
        (*rSvgFile) << "<circle class=\"node_" << region << "\" "
                    << "cx=\"" << scaled_x << "\" "
                    << "cy=\"" << scaled_y + mSvgSize << "\" "
                    << "r=\"" << rad << "\"/>" << "\n";
    }
    else if (scaled_y > mSvgSize - rad)
    {
        (*rSvgFile) << "<circle class=\"node_" << region << "\" "
                    << "cx=\"" << scaled_x << "\" "
                    << "cy=\"" << scaled_y - mSvgSize << "\" "
                    << "r=\"" << rad << "\"/>" << "\n";
    }
}

template <unsigned DIM>
void ThreeRegionSvgWriter<DIM>::AddGlyphToSvgFile(out_stream& rSvgFile,
                                                  c_vector<double, DIM> location,
                                                  unsigned region,
                                                  double rad,
                                                  double elongation,
                                                  int angle) const noexcept
{
    const double scaled_x = location[0] * mSvgSize;
    const double scaled_y = (1.0 - location[1]) * mSvgSize;

    const double max_size = std::max(rad, rad * elongation);

    (*rSvgFile) << "<ellipse class=\"glyph_" << region << "\" "
                << "transform=\"translate(" << scaled_x << " " << scaled_y
                << ") rotate(" << angle << ")\" "
                << "rx=\"" << rad << "\" "
                << "ry=\"" << elongation * rad << "\""
                << "/>" << "\n";

    // Account for possible wrap-around of glyph in x
    if (scaled_x < max_size)
    {
        (*rSvgFile) << "<ellipse class=\"glyph_" << region << "\" "
                    << "transform=\"translate(" << scaled_x + mSvgSize << " " << scaled_y
                    << ") rotate(" << angle << ")\" "
                    << "rx=\"" << rad << "\" "
                    << "ry=\"" << elongation * rad << "\""
                    << "/>" << "\n";
    }
    else if (scaled_x > mSvgSize - max_size)
    {
        (*rSvgFile) << "<ellipse class=\"glyph_" << region << "\" "
                    << "transform=\"translate(" << scaled_x - mSvgSize << " " << scaled_y
                    << ") rotate(" << angle << ")\" "
                    << "rx=\"" << rad << "\" "
                    << "ry=\"" << elongation * rad << "\""
                    << "/>" << "\n";
    }

    // Account for possible wrap-around of glyph in y
    if (scaled_y < max_size)
    {
        (*rSvgFile) << "<ellipse class=\"glyph_" << region << "\" "
                    << "transform=\"translate(" << scaled_x << " " << scaled_y + mSvgSize
                    << ") rotate(" << angle << ")\" "
                    << "rx=\"" << rad << "\" "
                    << "ry=\"" << elongation * rad << "\""
                    << "/>" << "\n";
    }
    else if (scaled_y > mSvgSize - max_size)
    {
        (*rSvgFile) << "<ellipse class=\"glyph_" << region << "\" "
                    << "transform=\"translate(" << scaled_x << " " << scaled_y - mSvgSize
                    << ") rotate(" << angle << ")\" "
                    << "rx=\"" << rad << "\" "
                    << "ry=\"" << elongation * rad << "\""
                    << "/>" << "\n";
    }
}

template <unsigned DIM>
void ThreeRegionSvgWriter<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template <unsigned DIM>
unsigned int ThreeRegionSvgWriter<DIM>::GetSamplingMultiple() const
{
    return mSamplingMultiple;
}

template <unsigned DIM>
void ThreeRegionSvgWriter<DIM>::SetSamplingMultiple(unsigned int samplingMultiple)
{
    mSamplingMultiple = samplingMultiple;
}

template <unsigned DIM>
double ThreeRegionSvgWriter<DIM>::GetSvgSize() const
{
    return mSvgSize;
}

template <unsigned DIM>
void ThreeRegionSvgWriter<DIM>::SetSvgSize(double svgSize)
{
    mSvgSize = svgSize;
}

template <unsigned DIM>
const std::array<unsigned int, 3>& ThreeRegionSvgWriter<DIM>::GetRegionSizes() const
{
    return mRegionSizes;
}

template <unsigned DIM>
void ThreeRegionSvgWriter<DIM>::SetRegionSizes(const std::array<unsigned int, 3>& regionSizes)
{
    mRegionSizes = regionSizes;
}

// Explicit instantiation
template class ThreeRegionSvgWriter<1>;
template class ThreeRegionSvgWriter<2>;
template class ThreeRegionSvgWriter<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ThreeRegionSvgWriter)
