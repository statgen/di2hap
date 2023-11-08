/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

#include <getopt.h>
#include <cstdlib>
#include <cmath>

std::vector<std::string> split_string_to_vector(const char* in, char delim)
{
  std::vector<std::string> ret;
  const char* d = nullptr;
  std::string token;
  const char* s = in;
  const char*const e = in + strlen(in);
  while ((d = std::find(s, e,  delim)) != e)
  {
    ret.emplace_back(std::string(s, d));
    s = d ? d + 1 : d;
  }
  ret.emplace_back(std::string(s,d));
  return ret;
}

class prog_args
{
private:
  std::vector<option> long_options_;
  std::string input_path_;
  std::string output_path_ = "/dev/stdout";
  std::string sex_map_path_;
  std::string haploid_code_ = "0";
  savvy::file::format output_format_ = savvy::file::format::sav;
  int compression_level_ = 6;
  bool verify_ = false;
  bool help_ = false;
  bool version_ = false;
public:
  prog_args() :
    long_options_(
      {
        {"haploid-code", required_argument, 0, 'c'},
        {"help", no_argument, 0, 'h'},
        {"output", required_argument, 0, 'o'},
        {"output-format", required_argument, 0, 'O'},
        {"sex-map", required_argument, 0, 'm'},
        {"version", no_argument, 0, 'v'},
        {"verify", no_argument, 0, 'V'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() const { return input_path_; }
  const std::string& output_path() const { return output_path_; }
  const std::string& sex_map_path() const { return sex_map_path_; }
  const std::string& haploid_code() const { return haploid_code_; }
  savvy::file::format output_format() const { return output_format_; }
  int compression_level() const { return compression_level_; }
  bool help_is_set() const { return help_; }
  bool version_is_set() const { return version_; }
  bool verify() const { return verify_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: di2hap [opts ...] input_file.{bcf,sav,vcf.gz} \n";
    os << "\n";
    os << " -c, --haploid-code   Code used for haploid samples in sex map (default: 0)\n";
    os << " -h, --help           Print usage\n";
    os << " -o, --output         Output path (default: /dev/stdout)\n";
    os << " -O, --output-format  Output file format (vcf, vcf.gz, bcf, ubcf, sav, usav; default: vcf)\n";
    os << " -m, --sex-map        Sex map file path (default: all samples are presumed haploid)\n";
    os << " -v, --version        Print version\n";
    os << " -V, --verify        Verify genotypes are homozygous before converting\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "c:hm:o:O:vV", long_options_.data(), &long_index)) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'c':
        haploid_code_ = optarg ? optarg : "";
        break;
      case 'h':
        help_ = true;
        return true;
      case 'o':
        output_path_ = optarg ? optarg : "";
        break;
      case 'O':
      {
        using fmt = savvy::file::format;
        std::string ot = optarg ? optarg : "";
        if (ot == "vcf")
        {
          output_format_ = fmt::vcf;
          compression_level_ = 0;
        }
        else if (ot == "vcf.gz")
        {
          output_format_ = fmt::vcf;
        }
        else if (ot == "bcf")
        {
          output_format_ = fmt::bcf;
        }
        else if (ot == "ubcf")
        {
          output_format_ = fmt::bcf;
          compression_level_ = 0;
        }
        else if (ot == "sav")
        {
          output_format_ = fmt::sav;
        }
        else if (ot == "usav")
        {
          output_format_ = fmt::sav;
          compression_level_ = 0;
        }
        else
        {
          std::cerr << "Invalid --output-format: " << ot << std::endl;
          return false;
        }
        break;
      }
      case 'm':
        sex_map_path_ = optarg ? optarg : "";
        break;
      case 'v':
        version_ = true;
        return true;
      case 'V':
        verify_ = true;
        break;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 0)
    {
      input_path_ = "/dev/stdin";
    }
    else if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
    }
    else
    {
      return std::cerr << "Error: invalid number of arguments\n", false;
    }

    return true;
  }
};

typedef std::int8_t gt_type;

bool verify(const std::vector<gt_type>& gt, const std::vector<int>& sex_map, const savvy::variant& rec, const std::vector<std::string>& sample_ids)
{
  std::size_t stride = gt.size() / sex_map.size();
  for (std::size_t i = 0; i < sex_map.size(); ++i)
  {
    if (!sex_map[i]) continue;

    for (std::size_t j = 1; j < stride; ++j)
    {
      if (gt[i * stride] != gt[i * stride + j])
      {
        std::cerr << "Error: cannot convert heterozygous to haploid at " << rec.chrom() << ":" << rec.pos() << ":" << rec.ref() << ":";
        for (auto it = rec.alts().begin(); it != rec.alts().end(); ++it)
        {
          if (it != rec.alts().begin())
            std::cerr << ",";
          std::cerr << *it;
        }
        std::cerr << ":" << sample_ids[i] << std::endl;
        return false;
      }
    }
  }

  return true;
}

int main(int argc, char** argv)
{
  prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }

  if (args.version_is_set())
  {
    std::cout << "hds-util v" << VERSION << std::endl;
    return EXIT_SUCCESS;
  }

  savvy::reader input_file(args.input_path());
  if (!input_file)
    return std::cerr << "Error: could not open input file\n", EXIT_FAILURE;

  savvy::writer output_file(args.output_path(), args.output_format(), input_file.headers(), input_file.samples(), args.compression_level());
  if (!output_file)
    return std::cerr << "Error: could not open output file\n", EXIT_FAILURE;

  std::vector<int> sex_map(input_file.samples().size(), 1);
  if (args.sex_map_path().size())
  {
    std::unordered_map<std::string, std::size_t> id_to_idx;
    id_to_idx.reserve(input_file.samples().size());
    for (std::size_t i = 0; i < input_file.samples().size(); ++i)
      id_to_idx[input_file.samples()[i]] = i;

    std::string line;
    std::ifstream sex_map_file(args.sex_map_path());
    while (std::getline(sex_map_file, line))
    {
      auto fields = split_string_to_vector(line.c_str(), '\t');
      if (fields.size() < 2)
        return std::cerr << "Error: malformed sex map\n", EXIT_FAILURE;

      auto res = id_to_idx.find(fields[0]);
      if (res == id_to_idx.end())
      {
        std::cerr << "Warning: Sex map ID not in VCF (" << fields[0] << ")" << std::endl;
      }
      else
      {
        if (fields[1] != args.haploid_code())
          sex_map[res->second] = 0;
      }
    }
  }

  int haploid_count = std::accumulate(sex_map.begin(), sex_map.end(), 0);
  std::cerr << "Notice: converting " << haploid_count << " samples to haploid" << std::endl;

  savvy::variant rec;
  std::vector<gt_type> gt;
  while (input_file >> rec)
  {
    rec.get_format("GT", gt);

    std::size_t stride = gt.size() / sex_map.size();

    if (haploid_count == input_file.samples().size())
    {
      if (args.verify() && !verify(gt, sex_map, rec, input_file.samples()))
        return EXIT_FAILURE;

      for (std::size_t i = 0; i < haploid_count; ++i)
        gt[i] = gt[i * stride];

      gt.resize(haploid_count);
    }
    else
    {
      if (args.verify() && !verify(gt, sex_map, rec, input_file.samples()))
        return EXIT_FAILURE;

      for (std::size_t i = 0; i < sex_map.size(); ++i)
      {
        if (sex_map[i])
        {
          for (std::size_t j = 1; j < stride; ++j)
            gt[i * stride + j] = savvy::typed_value::end_of_vector_value<gt_type>();
        }
      }
    }
   
    rec.set_format("GT", gt); 
    output_file << rec;
  }


  return input_file.bad() || !output_file.good() ? EXIT_FAILURE : EXIT_SUCCESS;
}


