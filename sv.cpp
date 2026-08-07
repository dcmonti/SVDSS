#include "sv.hpp"

using namespace std;

SV::SV() { l = 0; }

SV::SV(const string type_, const string &chrom_, uint s_, const string &refall_,
       const string &altall_, const uint w_, const uint cov_, const int ngaps_,
       const int score_, bool imprecise_, uint l_, string cigar_) {
  type = type_;
  chrom = chrom_;
  s = s_;
  refall = refall_;
  altall = altall_;
  // For symbolic alleles (<DEL>, <DUP>, <INV>) refall is just the anchor base;
  // END must be derived from l (SV length) so that END spans the affected region.
  // For sequence alleles (precise calls) refall contains the actual sequence.
  if (!altall_.empty() && altall_[0] == '<' &&
      (type_ == "DEL" || type_ == "DUP" || type_ == "INV") && l_ > 0) {
    e = s + (int)l_;
  } else {
    e = s + (int)refall.size() - 1;
  }
  w = w_;
  l = l_;
  cov = cov_;
  cov0 = 0;
  cov1 = 0;
  cov2 = 0;
  ngaps = ngaps_;
  score = score_;
  imprecise = imprecise_;
  cigar = cigar_;
  idx = type + "_" + chrom + ":" + to_string(s) + "-" + to_string(e);
  idx += "_" + to_string(abs(l));
  gt = "./.";
  rvec = "";
}

void SV::add_reads(const vector<string> &names) {
  for (const string &name : names) {
    if (name.empty())
      continue;
    if (!reads.empty())
      reads += ",";
    reads += name;
  }
}

void SV::add_sa_reads(const vector<string> &names) {
  for (const string &name : names) {
    if (name.empty())
      continue;
    if (!sa_reads.empty())
      sa_reads += ",";
    sa_reads += name;
  }
}

void SV::set_cov(int _cov, int _cov0, int _cov1, int _cov2) {
  cov = _cov;
  cov0 = _cov0;
  cov1 = _cov1;
  cov2 = _cov2;
}

void SV::set_rvec(const vector<tuple<int, int>> &reads) {
  for (const auto &r : reads)
    rvec += to_string(get<0>(r)) + ":" + to_string(get<1>(r)) + "-";
  rvec.pop_back();
}

void SV::set_gt(const string &_gt, int _gtq) {
  gt = _gt;
  gtq = _gtq;
}

ostream &operator<<(ostream &os, const SV &sv) {
  os << sv.chrom << "\t" << sv.s << "\t" << sv.idx << "\t" << sv.refall << "\t"
     << sv.altall << "\t"
     << "."
     << "\t"
     << "PASS"
     << "\t"
     // INFO
     << "VARTYPE=SV;"
     << "SVTYPE=" << sv.type << ";"
     << "SVLEN=" << (sv.type == "DEL" ? -sv.l : sv.l) << ";"
     << "END=" << sv.e << ";"
     << "WEIGHT=" << sv.w << ";"
     << "COV=" << sv.cov << ";"
     << "COV0=" << sv.cov0 << ";"
     << "COV1=" << sv.cov1 << ";"
     << "COV2=" << sv.cov2 << ";"
     << "AS=" << sv.score << ";"
     << "NV=" << sv.ngaps << ";"
     << "CIGAR=" << sv.cigar << ";";
  // BND breakends come in pairs; MATEID names the record holding the other end
  // so a consumer can collapse the two into one event without re-deriving the
  // junction from the ALT string. Emitted only when the mate was actually
  // written out.
  if (!sv.mate_id.empty())
    os << "MATEID=" << sv.mate_id << ";";
  // Junction detail, Manta/DRAGEN style. Emitted only when present, so a flush
  // junction produces exactly the record it did before these fields existed.
  if (sv.ins_len > 0) {
    os << "SVINSLEN=" << sv.ins_len << ";";
    if (!sv.ins_seq.empty())
      os << "SVINSSEQ=" << sv.ins_seq << ";";
  }
  if (sv.hom_len > 0) {
    os << "HOMLEN=" << sv.hom_len << ";";
    if (!sv.hom_seq.empty())
      os << "HOMSEQ=" << sv.hom_seq << ";";
  }
  os << "RVEC=" << sv.rvec << ";"
    << "READS=" << sv.reads << ";"
    << "SA_READS=" << sv.sa_reads
     << (sv.imprecise ? ";IMPRECISE\t" : "\t")
     // -
     << "GT:GQ"
     << "\t" << sv.gt << ":" << sv.gtq;
  return os;
}
