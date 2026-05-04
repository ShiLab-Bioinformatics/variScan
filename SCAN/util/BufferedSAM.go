package util


type SamBuffered struct {
        qf *Qfile
        currentReadName, currentChr1, currentCigar1, currentChr2, currentCigar2 string
        currentFlags1, currentPos1, currentFlags2, currentPos2 int
        currentExt1, currentExt2 []SAMAppendix
        isEOF bool
}

func (sb *SamBuffered) GetNextPair(){
        if sb.isEOF{
                sb.currentChr1 = ""
                sb.currentChr2 = ""
        }else{
                var err1 error
                var rn2 string

                sb.currentReadName, sb.currentFlags1, sb.currentChr1, sb.currentPos1, _, sb.currentCigar1, _, _, _, _, _, sb.currentExt1, err1 = sb.qf.SAMRecordEx()
                if err1 != nil{ sb.isEOF = true }
                rn2 , sb.currentFlags2, sb.currentChr2, sb.currentPos2, _, sb.currentCigar2, _, _, _, _, _, sb.currentExt2,   _  = sb.qf.SAMRecordEx()
                if rn2 != sb.currentReadName {panic("Unmatched names")}
        }
        return
}

func CreateBufferedSam(qs * Qfile) *SamBuffered{
        sb := new(SamBuffered)
        sb.qf = qs
        sb.GetNextPair()
        return sb
}

func (sb *SamBuffered) HasThisName(rname string) (hasThis bool, chr1, chr2 string, flag1, flag2, pos1, pos2 int, cigar1, cigar2 string, ex1, ex2 []SAMAppendix, isEOF bool) {
        if rname == sb.currentReadName{
                hasThis = true
                isEOF = sb.isEOF

                chr1 = sb.currentChr1
                chr2 = sb.currentChr2
                flag1 = sb.currentFlags1
                flag2 = sb.currentFlags2
                cigar1 = sb.currentCigar1
                cigar2 = sb.currentCigar2
                ex1 = sb.currentExt1
                ex2 = sb.currentExt2
                pos1 = sb.currentPos1
                pos2 = sb.currentPos2

                sb.GetNextPair()
                return
        }else{
                hasThis = false
                return
        }
}

