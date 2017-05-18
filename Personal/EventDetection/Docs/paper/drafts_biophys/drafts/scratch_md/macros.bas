Attribute VB_Name = "NewMacros"
Sub wrap_table()
Attribute wrap_table.VB_ProcData.VB_Invoke_Func = "Normal.NewMacros.Macro3"
'
' wrap_table Macro
'
'
    With ActiveDocument.Tables(1).Rows
        .WrapAroundText = True
        .HorizontalPosition = wdTableLeft
        .RelativeHorizontalPosition = wdRelativeHorizontalPositionColumn
        .DistanceLeft = InchesToPoints(0.13)
        .DistanceRight = InchesToPoints(0.13)
        .VerticalPosition = InchesToPoints(2.1)
        .RelativeVerticalPosition = wdRelativeVerticalPositionPage
        .DistanceTop = InchesToPoints(0)
        .DistanceBottom = InchesToPoints(0)
        .AllowOverlap = False
    End With
    Selection.ParagraphFormat.Alignment = wdAlignParagraphJustify
End Sub
